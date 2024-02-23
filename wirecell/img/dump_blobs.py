#!/usr/bin/env python3
'''
Dump out signatures of blobs for debugging
signature: [tmin, tmax, umin, umax, vmin, vmax, wmin, wmax]
'''
from wirecell import units
import matplotlib.pyplot as plt
import numpy
import math
def quadrature_sum(lst):
    squares_sum = sum(x**2 for x in lst)
    return math.sqrt(squares_sum)

def bsignature(gr, bnode, tick=500, focus='val'):
    '''
    focus: val or unc
    '''
    sig = []
    id2name = {1:'u', 2:'v', 4:'w'}
    chan_index = dict()
    chan_status = dict()
    signal = dict()
    meas = dict()
    for id in id2name:
        chan_index[id] = []
        chan_status[id] = []
    for node in gr.neighbors(bnode):
        ndata = gr.nodes[node]
        if ndata['code'] == 'm':
            # print(ndata)
            wpid = ndata['wpid']
            val = ndata['val']
            meas[wpid] = val
    for node in gr.neighbors(bnode):
        ndata = gr.nodes[node]
        if ndata['code'] == 's':
            # print(ndata)
            tmin = ndata['start']//tick
            tmax = tmin + ndata['span']//tick
            sig.append(tmin)
            sig.append(tmax)
            signal.update({s['ident']:s for s in ndata['signal']})

    for node in gr.neighbors(bnode):
        ndata = gr.nodes[node]
        if ndata['code'] == 'w':
            # print(ndata)
            # chid: global; index: per-plane
            chid = ndata['chid']
            wpid = ndata['wpid']
            index = ndata['index']
            chan_index[wpid].append(index)
            if chid in signal:
                val = signal[chid][focus]
            else:
                val = -1
                # for key in sorted(signal):
                #     print(key, ': ', signal[key]['val'])
                # raise RuntimeError(f'{chid} not in signal')
            if val < 1:
                val = 0
            chan_status[wpid].append(val)
    chan_offset = {1:0, 2:2400, 4: 4800}
    for wpid in chan_index:
        # print(wpid, chan_index[wpid])
        if len(chan_index[wpid]) == 0:
            # if len(sig) == 0 or sig[0] == 1024:
            #     print(gr.nodes[bnode])
            #     for node in gr.neighbors(bnode):
            #         ndata = gr.nodes[node]
            #         print(ndata)
            #     print('')
            return None
        min = numpy.min(chan_index[wpid]) + chan_offset[wpid]
        max = numpy.max(chan_index[wpid]) + chan_offset[wpid]
        sig.append(min)
        sig.append(max)
    for wpid in chan_status:
        if focus not in ['val', 'unc']:
            raise('focus not in [\'val\', \'unc\']')
        if focus == 'val':
            sig.append(int(sum(chan_status[wpid])))
        if focus == 'unc':
            sig.append(int(quadrature_sum(chan_status[wpid])))
    for wpid in id2name:
        if wpid in meas:
            sig.append(int(meas[wpid]))
        else:
            sig.append(int(0))
    # sig.append(gr.nodes[bnode]['value'])
    # sig.append(gr.nodes[bnode]['ident'])
    # print(f'sig: {sig}')
    return sig

def _sort(arr):
    ind = numpy.lexsort((arr[:,7],arr[:,6],arr[:,5],arr[:,4],arr[:,3],arr[:,2],arr[:,1],arr[:,0]))
    arr = numpy.array([arr[i] for i in ind])
    return arr

def dump_blobs(gr, bvalscale=1.0/0.8, sigfile=None, dumpfile="/dev/stdout"):
    '''
    extract blob signatures
    '''
    out = open(dumpfile, "w")
    sigs = []
    count = 0
    for node, ndata in gr.nodes.data():
        if ndata['code'] != 'b':
            continue;
        sig = bsignature(gr, node)
        bval = ndata['val'] * bvalscale
        sig.append(int(bval))
        sig.append(ndata['ident'])
        if sig is not None:
            sigs.append(sig)
        else:
            count += 1
    sigs = numpy.array(sigs)
    print('WCT:')
    print(f'sigs.shape: {sigs.shape}')
    # sigs = sigs[sigs[:,8]>0,:]
    # sigs = sigs[sigs[:,9]>0,:]
    # sigs = sigs[sigs[:,10]>0,:]
    sigs = _sort(sigs)
    # for i in range(min([sigs.shape[0], 20])):
    numpy.save("wct.npy", sigs[:,:15])
    for i in range(sigs.shape[0]):
        # print(i, sigs[i,:])
        print(sigs[i,0:2],                    # tick
            sigs[i,2], ':', sigs[i,3]+1, ',', # u wire bounds
            sigs[i,4], ':', sigs[i,5]+1, ','  # v wire bounds
            ,sigs[i,6], ':', sigs[i,7]+1      # w wire bounds
            ,sigs[i,8:11]                     # sum of wire charge
            ,sigs[i,11:14]                    # measurement
            ,sigs[i,14]                       # blob charge
            # ,sigs[i,15]                       # blob ident
            )

    print(f'sum > -1 {sum(sigs[sigs[:,14]>-1,14])}')
    print(f'sum > 300 {sum(sigs[sigs[:,14]>300,14])}')

    print(f'> -1 {len(sigs[sigs[:,14]>-1,:])}')
    print(f'> 300 {len(sigs[sigs[:,14]>300,:])}')
    print(f'> 1000 {len(sigs[sigs[:,14]>1000,:])}')
    print(f'> 10000 {len(sigs[sigs[:,14]>10000,:])}')
    print(f'> 100000 {len(sigs[sigs[:,14]>100000,:])}')
    if sigfile is not None:
        numpy.save(sigfile, sigs)