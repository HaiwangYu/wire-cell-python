#!/usr/bin/python3
'''
Fixme: make this into a proper click main
'''
import os
import sys
import json
from collections import defaultdict
import subprocess

import click

from wirecell import units
from wirecell.util import jsio
from wirecell.util.cli import jsonnet_loader

cmddef = dict(context_settings = dict(help_option_names=['-h', '--help']))

@click.group("pgraph", **cmddef)
@click.pass_context
def cli(ctx):
    '''
    Wire Cell Signal Processing Features
    '''

class Node (object):
    def __init__(self, tn, params=True, **attrs):
        if not attrs:
            print ("Node(%s) with no attributes"%tn)

        self.tn = tn
        self._params = params
        tn = tn.split(":")
        self.type = tn[0]
        try:
            self.name = tn[1]
        except IndexError:
            self.name = ""
        self.ports = defaultdict(set);
        self.attrs = attrs

    @property
    def display_name(self):
        if self.name:
            return "[%s]"%self.name
        return "(unnamed)"

    def add_port(self, end, ident):
        '''
        Add a port of end "head" or "tail" and ident (number).
        '''
        self.ports[end].add(ident)

    def dot_name(self, port=None):
        return self.tn.replace(":","_")

    def dot_label_one(self, v, recur=True):
        if isinstance(v,list):
            siz = len(v)
            psize = min(siz, 3)
            last = ""
            if siz > psize:
                last = "..."
            if recur:
                vstr = ",".join([self.dot_label_one(vv, False) for vv in v[:psize]])
            if not recur or psize < siz:
                vstr = ",..."
            v = "list(%d):[%s]"%(siz, vstr)
            return v
        if isinstance(v,dict):
            v = "dict(%d):[%s]"%(len(v), self.dot_label_one(list(v.keys()), False))
            return v
        return str(v)


    def dot_label(self):
        ret = list()
        if "head" in self.ports:
            head = "{%s}" % ("|".join(["<in%d>%d"%(num,num) for num in sorted(self.ports["head"])]),)
            ret.append(head)

        body = [self.type, self.display_name]
        if self._params:
            for k,v in sorted(self.attrs.items()):
                v = self.dot_label_one(v)
                one = "%s = %s" % (k,v)
                body.append(one)
        body = r"\n".join(body)
        body = r"{%s}" % body
        ret.append(body)

        if "tail" in self.ports:
            tail = "{%s}" % ("|".join(["<out%d>%d"%(num,num) for num in sorted(self.ports["tail"])]),)
            ret.append(tail)

        return "{%s}" % ("|".join(ret),)


def is_string(x):
    return type(x) in [type(u""), type("")]
def is_list(x):
    return type(x) in [list]
def is_list_of_string(x):
    if not is_list(x): return False
    return all(map(is_string, x))

def dotify(edge_dat, attrs, params=True, services=True):
    '''
    Return GraphViz text.

    If attrs is a dictionary, append to the node a list of its items.

    If params is True, show the attributes.

    If services is True, include non DFP node components.

    '''

    nodes = dict()
    def get(edge, end):
        tn = edge[end]["node"]
        try:
            n = nodes[tn]
        except KeyError:
            n = Node(tn, params, **attrs.get(tn, {}))
            nodes[tn] = n
        p = edge[end].get("port",0)
        n.add_port(end, p)
        return n,p
    
    edges = list()
    for edge in edge_dat:
        t,tp = get(edge, "tail")
        h,hp = get(edge, "head")
        e = '"%s":out%d -> "%s":in%d' % (t.dot_name(),tp, h.dot_name(),hp)
        edges.append(e);

    # Try to find non DFP node components referenced.
    if services:
        for tn,n in list(nodes.items()):
            for k,v in n.attrs.items():
                tocheck = None
                if is_string(v):
                    tocheck = [v]
                if is_list_of_string(v):
                    tocheck = v
                if not tocheck:
                    continue
                for maybe in tocheck:
                    if maybe not in attrs:
                        continue

                    cn = nodes.get(maybe,None);
                    if cn is None:
                        cn = Node(maybe, params, **attrs.get(maybe, {}))
                        nodes[maybe] = cn

                    e = '"%s" -> "%s"[style=dashed,color=gray]' % (n.dot_name(), cn.dot_name())
                    edges.append(e)

    ret = ["digraph pgraph {",
           "rankdir=LR;",
           "\tnode[shape=record];"]
    for nn,node in sorted(nodes.items()):
        nodestr = '\t"%s"[label="%s"];' % (node.dot_name(), node.dot_label())
        #print(nodestr)
        ret.append(nodestr)
    for e in edges:
        ret.append("\t%s;" % e)
    ret.append("}")
    return '\n'.join(ret);


# def jsonnet_try_path(path, rel):
#     if not rel:
#         raise RuntimeError('Got invalid filename (empty string).')
#     if rel[0] == '/':
#         full_path = rel
#     else:
#         full_path = os.path.join(path, rel)
#     if full_path[-1] == '/':
#         raise RuntimeError('Attempted to import a directory')

#     if not os.path.isfile(full_path):
#         return full_path, None
#     with open(full_path) as f:
#         return full_path, f.read()


# def jsonnet_import_callback(path, rel):
#     paths = [path] + os.environ.get("WIRECELL_PATH","").split(":")
#     for maybe in paths:
#         try:
#             full_path, content = jsonnet_try_path(maybe, rel)
#         except RuntimeError:
#             continue
#         if content:
#             return full_path, content
#     raise RuntimeError('File not found')



def resolve_path(obj, jpath):
    '''
    Select out a part of obj based on a "."-separated path.  Any
    element of the path that looks like an integer will be cast to
    one assuming it indexes an array.
    '''
    jpath = jpath.split('.')
    for one in jpath:
        if not one:
            break
        try:
            one = int(one)
        except ValueError:
            pass
        obj = obj[one]

    return obj

def uses_to_params(uses):
    '''
    Given a list of nodes, return a dictionary of their "data" entries
    keyed by 'type' or 'type:name'
    '''
    ret = dict()
    for one in uses:
        if type(one) != dict:
            print (type(one),one)
        tn = one[u"type"]
        if "name" in one and one['name']:
            # print (one["name"])
            tn += ":" + str(one["name"])
        ret[tn] = one.get("data", {})
    return ret

# fixme: add tla related options so any .jsonnet can also be loaded
@cli.command("dotify")
@click.option("--dpath", default="-1",
              help="A dot-delimited path into the data structure to locate a graph-like object")
@click.option("--params/--no-params", default=True,
              help="Enable/disable the inclusion of contents of configuration parameters") 
@click.option("--services/--no-services", default=True,
              help="Enable/disable the inclusion 'service' (non-node) type components") 
@jsonnet_loader("in-file")
@click.argument("out-file")
@click.pass_context
def cmd_dotify(ctx, dpath, params, services, in_file, out_file):
    '''Convert a WCT cfg to a GraphViz dot or rendered file.

    The config file may be JSON of Jsonnet.  It is expected to be in
    the form of a WCT configuration sequence with one (typically last)
    entry holding an object with a .data.edges key (eg a TbbFlow or
    Pgraph config object).  If a non-standard structure is given, the
    object with .edges can be located with --dpath

    Example bash command assuming WIRECELL_PATH properly set

    $ wirecell-pgraph dotify mycfg.jsonnet mycfg.pdf

    Or piecewise

    $ wcsonnet mycfg.jsonnet > mycfg.json

    $ wirecell-pgraph dotify mycfg.json mycfg.dot

    $ dot -Tpdf -o mycfg.pdf mycfg.dot

    If jsonnet is given the usual: -A/--tla, -J/--jpath args can be
    given.
    '''
    dat = in_file
    try: 
        cfg = resolve_path(dat, dpath)
    except Exception:
        click.echo('failed to resolve path "%s" in object:\n' % (dpath))
        sys.exit(1)

    # if cfg["type"] not in ["Pgrapher", "Pnode"]:
    #     click.echo('Object must be of "type" Pgrapher or Pnode, got "%s"' % cfg["type"])
    #     sys.exit(1)

    graph_apps = ("Pgrapher", "TbbFlow")

    # the Pgrapher app holds edges in "data" attribute
    if cfg.get("type","") in graph_apps:
        print ('Pgrapher object found at dpath: "%s" with %d nodes' % (dpath, len(dat)))
        edges = cfg["data"]["edges"] 
        uses = dat # if Pgrapher, then original is likely the full config sequence.
    else:
        edges = cfg["edges"] # Pnodes have edges as top-level attribute
        uses = cfg.get("uses", list())

    attrs = uses_to_params(uses)
    dtext = dotify(edges, attrs, params, services)
    ext = os.path.splitext(out_file)[1][1:]
    dot = "dot -T %s -o %s" % (ext, out_file)
    proc = subprocess.Popen(dot, shell=True, stdin = subprocess.PIPE)
    proc.communicate(input=dtext.encode("utf-8"))
    return

def main():
    cli(obj=dict())

if '__main__' == __name__:
    main()
    
