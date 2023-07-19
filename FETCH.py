#e.g. $python FETCH.py -f example -p my_cool_project -s FC114_A2_A02_002.fcs FC114_C1_C01_025.fcs

import argparse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import scipy

def _centered(arr, newsize):
    # Return the center newsize portion of the array.
    newsize = np.asarray(newsize)
    currsize = np.array(arr.shape)
    startind = (currsize - newsize) // 2
    endind = startind + newsize
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]


scipy.signal.signaltools._centered = _centered

import flowkit as fk
import numpy as np
from numpy.linalg import norm
import scipy.stats as st
from os.path import join
import os
#these are still libraries

from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from numpy.linalg import eig, inv
from matplotlib.patches import Ellipse
from matplotlib.path import Path
from sklearn.model_selection import GridSearchCV
from matplotlib.offsetbox import AnchoredText
from matplotlib import path
import pandas as pd 
import seaborn as sns

plt.switch_backend('agg')
plt.ioff()
#plot formatting

#def = define function; in this case... allows argument input into terminal.. 
def parse_arguments():
    ap = argparse.ArgumentParser(description="Parse arguments")

    ap.add_argument(
        "-f", "--folder",
        default="",
        help="The path to a directory containing your .fcs files",
        type=str )

    ap.add_argument(
        "-p", "--project",
        default="Untitled_project",
        help="The project name",
        type=str )

    ap.add_argument(
        "-s", "--skip_renaming",
        default='',
        help="This script expects fcs files to be formatted as 'FC114_A1_A01_001.fcs'; if just a few files are not formatted like that, pass them as arguments in this function to avoid errors; otherwise edit summarize function too parse your filenames correctly", nargs='*')

    return ap.parse_args()

#this function writes gatefiles that can be opened in flowjo specifically
def gate_writer(vertices1, vertices2, boundaries, filename, channame1=None, channame2=None, fluorophore1=None, fluorophore2=None):
    gate_text = ['<?xml version="1.0" encoding="UTF-8"?>',
    '<gating:Gating-ML',
    '    xmlns:gating="http://www.isac-net.org/std/Gating-ML/v2.0/gating"',
    '    xmlns:data-type="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes">']
    polygon1 = ['    <gating:PolygonGate gating:id="Polygon1">',
    '        <gating:dimension gating:compensation-ref="uncompensated">',
    '            <data-type:fcs-dimension data-type:name="SSC-A" />',
    '        </gating:dimension>',
    '        <gating:dimension gating:compensation-ref="uncompensated">',
    '            <data-type:fcs-dimension data-type:name="FSC-A" />',
    '        </gating:dimension>']
    for vertex in vertices1:
        vertex = ['        <gating:vertex>',
        '            <gating:coordinate data-type:value="' + str(vertex[0]) + '" />',
        '            <gating:coordinate data-type:value="' + str(vertex[1]) + '" />',
        '        </gating:vertex>']
        polygon1 = polygon1 + vertex
    polygon1 = polygon1 + ['    </gating:PolygonGate>']    
    if vertices2 is not None:
        polygon2 = ['    <gating:PolygonGate gating:id="Polygon2">',
        '        <gating:dimension gating:compensation-ref="uncompensated">',
        '            <data-type:fcs-dimension data-type:name="FSC-A" />',
        '        </gating:dimension>',
        '        <gating:dimension gating:compensation-ref="uncompensated">',
        '            <data-type:fcs-dimension data-type:name="FSC-H" />',
        '        </gating:dimension>']
        for vertex in vertices2:
            vertex = ['        <gating:vertex>',
            '            <gating:coordinate data-type:value="' + str(vertex[0]) + '" />',
            '            <gating:coordinate data-type:value="' + str(vertex[1]) + '" />',
            '        </gating:vertex>']
            polygon2 = polygon2 + vertex  
        polygon2 = polygon2 + ['    </gating:PolygonGate>']  
        and_gate = ['    <gating:BooleanGate gating:id="And1">',
        '        <data-type:custom_info>',
        '            Only keep results satisfying both Polygon gates',
        '        </data-type:custom_info>',
        '        <gating:and>',
        '            <gating:gateReference gating:ref="Polygon1" />',
        '            <gating:gateReference gating:ref="Polygon2" />',
        '        </gating:and>',
        '    </gating:BooleanGate>']
        if boundaries is not None:
            quadrants = ['    <gating:QuadrantGate gating:id="Quadrant1" gating:parent_id="And1">',
            '        <gating:divider gating:id="A" gating:compensation-ref="uncompensated">',
            '            <data-type:fcs-dimension data-type:name="' + channame1 + '" />',
            '            <gating:value>' + str(boundaries[1]) + '</gating:value>',
            '        </gating:divider>',
            '        <gating:divider gating:id="B" gating:compensation-ref="uncompensated">',
            '            <data-type:fcs-dimension data-type:name="' + channame2 + '" />',
            '            <gating:value>' + str(boundaries[0]) + '</gating:value>',
            '        </gating:divider>',
            '        <gating:Quadrant gating:id="Untransfected">',
            '            <gating:position gating:divider_ref="A" gating:location="' + str(boundaries[1] - 1) + '" />',
            '            <gating:position gating:divider_ref="B" gating:location="' + str(boundaries[0] - 1) + '" />',
            '        </gating:Quadrant>',
            '        <gating:Quadrant gating:id="Double-Positive">',
            '            <gating:position gating:divider_ref="A" gating:location="' + str(boundaries[1] + 1) + '" />',
            '            <gating:position gating:divider_ref="B" gating:location="' + str(boundaries[0] + 1) + '" />',
            '        </gating:Quadrant>',
            '        <gating:Quadrant gating:id="' + fluorophore1 + '">',
            '            <gating:position gating:divider_ref="A" gating:location="' + str(boundaries[1] + 1) + '" />',
            '            <gating:position gating:divider_ref="B" gating:location="' + str(boundaries[0] - 1) + '" />',
            '        </gating:Quadrant>',
            '        <gating:Quadrant gating:id="' + fluorophore2 + '">',
            '            <gating:position gating:divider_ref="A" gating:location="' + str(boundaries[1] - 1) + '" />',
            '            <gating:position gating:divider_ref="B" gating:location="' + str(boundaries[0] + 1) + '" />',
            '        </gating:Quadrant>',
            '    </gating:QuadrantGate>']
            second_and_gate = ['    <gating:BooleanGate gating:id="And2Untransfected">',
            '        <data-type:custom_info>',
            '            Only keep results satisfying both Polygon gates and Untransfected',
            '        </data-type:custom_info>',
            '        <gating:and>',
            '            <gating:gateReference gating:ref="And1" />',
            '            <gating:gateReference gating:ref="Untransfected" />',
            '        </gating:and>',
            '    </gating:BooleanGate>']
            third_and_gate = ['    <gating:BooleanGate gating:id="And3DoublePositive">',
            '        <data-type:custom_info>',
            '            Only keep results satisfying both Polygon gates and Double-Positive',
            '        </data-type:custom_info>',
            '        <gating:and>',
            '            <gating:gateReference gating:ref="And1" />',
            '            <gating:gateReference gating:ref="Double-Positive" />',
            '        </gating:and>',
            '    </gating:BooleanGate>']
            fourth_and_gate = ['    <gating:BooleanGate gating:id="And4' + fluorophore1 + '">',
            '        <data-type:custom_info>',
            '            Only keep results satisfying both Polygon gates and ' + fluorophore1,
            '        </data-type:custom_info>',
            '        <gating:and>',
            '            <gating:gateReference gating:ref="And1" />',
            '            <gating:gateReference gating:ref="' + fluorophore1 + '" />',
            '        </gating:and>',
            '    </gating:BooleanGate>']
            fifth_and_gate = ['    <gating:BooleanGate gating:id="And5' + fluorophore2 + '">',
            '        <data-type:custom_info>',
            '            Only keep results satisfying both Polygon gates and ' + fluorophore2,
            '        </data-type:custom_info>',
            '        <gating:and>',
            '            <gating:gateReference gating:ref="And1" />',
            '            <gating:gateReference gating:ref="' + fluorophore2 + '" />',
            '        </gating:and>',
            '    </gating:BooleanGate>']
            gate_text = gate_text + polygon1 + polygon2 + and_gate + quadrants + \
            second_and_gate + third_and_gate + fourth_and_gate + fifth_and_gate + ['</gating:Gating-ML>']
        else:
            gate_text = gate_text + polygon1 + polygon2 + and_gate + ['</gating:Gating-ML>']
    else: 
        gate_text = gate_text + polygon1 + ['</gating:Gating-ML>']
    with open(filename, "w") as f:
        for line in gate_text:
            if line not in ['\n', '\r\n']:
                f.write("%s\n" % line)

#selection for the first gate
def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    b,c,d,f,g,a=a[1]/2., a[2], a[3]/2., a[4]/2., a[5], a[0]
    num=b*b-a*c
    cx=(c*d-b*f)/num
    cy=(a*f-b*d)/num
    
    angle=0.5*np.arctan(2*b/(a-c))*180/np.pi
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    a=np.sqrt(abs(up/down1))
    b=np.sqrt(abs(up/down2))
    ell=Ellipse((cx,cy),a*2.,b*2.,angle=angle)
    ell_coord=ell.get_verts()
    return ell_coord

#kde = kernel density estimation; matching density to color    
def make_kde(points):
    x = points[:, 0]
    y = points[:, 1]
    # Define the borders
    deltaX = (max(x) - min(x))/1000
    deltaY = (max(y) - min(y))/1000
    xmin = min(x) - deltaX
    xmax = max(x) + deltaX
    ymin = min(y) - deltaY
    ymax = max(y) + deltaY
    # Create meshgrid
    xx, yy = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    fig = plt.figure(figsize=(16,16))
    ax = fig.gca()
    cset1 = ax.contour(xx, yy, f, levels=50, colors='k')
    figure_centre = [(xmin + xmax)/2, (ymin + ymax)/2]
    plt.cla()
    plt.clf()
    return [cset1.allsegs, figure_centre, x, y]


#similar to above
def getKernelDensityEstimation(values, x, bandwidth = 0.2, kernel = 'gaussian'):
    model = KernelDensity(kernel = kernel, bandwidth=bandwidth)
    model.fit(values[:, np.newaxis])
    log_density = model.score_samples(x[:, np.newaxis])
    return np.exp(log_density)

#a helper function for kde color matching
def bestBandwidth(data, minBandwidth = 0.1, maxBandwidth = 2, nb_bandwidths = 30, cv = 30):
    """
    Run a cross validation grid search to identify the optimal bandwidth for the kernel density
    estimation.
    """
    model = GridSearchCV(KernelDensity(),
                        {'bandwidth': np.linspace(minBandwidth, maxBandwidth, nb_bandwidths)}, cv=cv, n_jobs=-1) 
    model.fit(data)
    return model.best_params_['bandwidth']


#possibly for defining the last gate
def z(samplename, Z, a, vertices1, vertices2, dest, sample, fluorophore1, fluorophore2, channame1, channame2):
    fluorophores = fluorophore1 + '_' + fluorophore2
    filename = join(dest, 'gates.xml')
    d = Z[a]
    if len(d) == 0:
        return [samplename, 0, None, 0]
    new_Z = np.array([Z[a][:, 0] + abs(min(Z[a][:, 0])) + 1, Z[a][:, 1] + abs(min(Z[a][:, 1])) + 1])
    new_Z= new_Z.T
    Z_log = np.log(new_Z) 
    try:
        [alls, figure_centre, x, y] = make_kde(Z_log)
    except ValueError:
        return [samplename, 0, None, 0]
    max_area = 0
    best_top_point = None
    best_right_point = None
    
    xy = np.vstack([x,y])
    cv_bandwidth = bestBandwidth(xy.T)
    kde_model = KernelDensity(kernel='gaussian', bandwidth=cv_bandwidth).fit(xy.T)
    kde = np.exp(kde_model.score_samples(xy.T))
    idx = kde.argsort()
    
    candidates = []
    for j in range(len(alls)):
        for ii, seg in enumerate(alls[j]):
            p = Path(seg) # make a polygon
            grid = p.contains_points(xy.T)
            mean_kde = np.mean(kde[grid])
            top_point_log = max(seg[:,1])
            right_point_log = max(seg[:,0])
            top_point = np.exp(top_point_log) - abs(min(Z[a][:, 1])) - 1
            right_point = np.exp(right_point_log) - abs(min(Z[a][:, 0])) - 1
            transfected_cells_x = right_point >= 500
            transfected_cells_y = top_point >= 500
            plt.plot(seg[:,0], seg[:,1], '.-')
            if transfected_cells_x or transfected_cells_y:
                continue 
            area = 0.5*np.abs(np.dot(seg[:,0],np.roll(seg[:,1],1))-np.dot(seg[:,1],np.roll(seg[:,0],1)))    
            if area > max_area:
                candidates.append([area, kde[grid], top_point, right_point])
    largest_cand = [0]
    if len(candidates) == 0:
        print("Pipeline error on this file")
        return [samplename, 0, None, 0]
    for cand in candidates:
        if cand[0] > largest_cand[0]:
            largest_cand = cand
    h_plt = plt.hist(largest_cand[1], bins=100)
    bin_count = h_plt[0]
    cutoff = h_plt[1]        
    quantile_cutoff = np.quantile(cutoff, 0.60)
    
    for cand in candidates:
        if np.mean(cand[1]) < quantile_cutoff:
                continue
        if cand[0] > max_area:
                max_area = cand[0]
                best_top_point = cand[2]
                best_right_point = cand[3]
      
    
    boundaries = [best_right_point, best_top_point]
    if len(sample.channels) == 7:
        gname = join(dest, fluorophores + '_gates.xml')
    else:
        gname = join(dest, 'gates.xml')
    gate_writer(vertices1, vertices2, boundaries, gname, channame1, channame2, fluorophore1, fluorophore2)
    g_strat = fk.parse_gating_xml(gname)    
    gs_results = g_strat.gate_sample(sample)
    e = gs_results.get_gate_membership('And3DoublePositive')
    fig = figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k') 
    
    xy = d
    x = d[:, 0]
    y = d[:, 1]
    x_sorted, y_sorted, kde_sorted = x[idx], y[idx], kde[idx]
    # parameters of the main output plot
    plt.scatter(x, y, c=kde, cmap = 'turbo', s=15)
    plt.yscale('symlog', linthresh=1000)
    plt.xscale('symlog', linthresh=1000)
    #these are gate lines; min, max are the range of point values; best points define position of the gate
    plt.plot([min(x), max(x)], [best_top_point, best_top_point], c='black')
    plt.plot([best_right_point, best_right_point], [min(y), max(y)], c='black')
    #df is a table format for... parsed.. data..
    df = gs_results.report
    df = df.reset_index()
    #numbers of each individual quadrant
    double_positives = list(df.loc[df['gate_name'] == 'And3DoublePositive']['count'])[0]
    green = list(df.loc[df['gate_name'] == 'And5' + fluorophore2]['count'])[0]
    red = list(df.loc[df['gate_name'] == 'And4' + fluorophore1]['count'])[0]
    untransfected = list(df.loc[df['gate_name'] == 'And2Untransfected']['count'])[0]
    #this is a contingency for blank samples
    if double_positives + green + red == 0:
        return [samplename, 0, None, 0]
    #defining the FETCH  score
    FETCH_score = double_positives/(double_positives + green + red)
    #another contingency
    if FETCH_score > 0.90 or untransfected/(double_positives + green + red + untransfected) > 0.90:
        return [samplename, 0, None, double_positives + green + red + untransfected]
    r_g = red/green
    ax = plt.gca()
    minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    ax.yaxis.set_minor_locator(minor)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.xaxis.set_minor_locator(minor)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(which='minor', length=10, width=2)
    ax.tick_params(which='major', length=20, width=3)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=14)
    plt.setp(ax.get_yticklabels(), fontsize=14)
    #text boxes in the output plot
    txt1 = AnchoredText('Q1\n' + str(round(100*red/(double_positives + green + red + untransfected), 1)), loc="upper left", pad=0.4, borderpad=0, prop={"fontsize":14})
    txt2 = AnchoredText('Q2\n' + str(round(100*double_positives/(double_positives + green + red + untransfected), 1)), loc="upper right", pad=0.4, borderpad=0, prop={"fontsize":14})
    txt3 = AnchoredText('Q3\n' + str(round(100*green/(double_positives + green + red + untransfected), 1)), loc="lower right", pad=0.4, borderpad=0, prop={"fontsize":14})
    txt4 = AnchoredText('Q4\n' + str(round(100*untransfected/(double_positives + green + red + untransfected), 1)), loc="lower left", pad=0.4, borderpad=0, prop={"fontsize":14})
    #this puts texts boxes onto the plot
    ax.add_artist(txt1)
    ax.add_artist(txt2)
    ax.add_artist(txt3)
    ax.add_artist(txt4)
    ax.set_title("FETCH Score: " + str(FETCH_score))
    #look into the bbox, bounding box
    plt.savefig(join(dest, fluorophores + '_double_positive_final.pdf'), format='pdf', bbox_inches='tight')
    #clear axes and figure to plot next
    plt.cla()
    plt.clf()
    return [samplename, FETCH_score, r_g, double_positives + green + red + untransfected]

 #Above, the helper files were added, Below here, the actual processing, central functions are listed
def FETCH_analysis(inputlist):
    plt.close('all')
    fcs_path, samplename, dest = inputlist
    #make a directly for an FCS file and use flowkit to parse that, to get variable called sample
    os.mkdir(dest)
    plt.grid(visible=None)
    sample = fk.Sample(fcs_path)
    arr1 = sample.get_channel_events(0, source='raw', subsample=False) #FSC-A
    arr2 = sample.get_channel_events(2, source='raw', subsample=False) #SSC-A
    arr3 = sample.get_channel_events(1, source='raw', subsample=False) #FSC-H
    arr4 = sample.get_channel_events(3, source='raw', subsample=False) #1-A or FITC-A(Emerald)
    if len(sample.channels) == 6:
        arr5 = sample.get_channel_events(4, source='raw', subsample=False) #5-A(RFP670)
        Z = np.stack((arr4, arr5), axis=1)
    elif len(sample.channels) == 5:
        arr5 = sample.get_channel_events(4, source='raw', subsample=False) #PE-Texas Red-A(mCherry)
        Z = np.stack((arr4, arr5), axis=1)
    elif len(sample.channels) == 7:
        arr5 = sample.get_channel_events(4, source='raw', subsample=False) #PE-A (mApple)
        arr6 = sample.get_channel_events(5, source='raw', subsample=False) #APC-A (RFP670)
        Z_ea = np.stack((arr4, arr5), axis=1)
        Z_er = np.stack((arr4, arr6), axis=1)
        Z_ar = np.stack((arr5, arr6), axis=1)
    else:
        raise Exception("Something is wrong with your channel number")
        #arr = array, plot.. x = 1st gate and y = 2nd gate
    X = np.stack((arr2, arr1), axis=1)
    Y = np.stack((arr1, arr3), axis=1)     
    #this loop goes through the contours of the first gate 
    [alls, figure_centre, x, y] = make_kde(X)
    max_area = 0
    best_seg = None
    point_num = None
    plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k')
    seg_list = []
    min_points = 4000
    for j in range(len(alls)):
        for ii, seg in enumerate(alls[j]):
            non_single_cells_x = min(seg[:,0]) <= 25000
            non_single_cells_y = min(seg[:,1]) <= 25000
            out_of_bounds_x = max(seg[:,0]) > (max(x) - 10000)
            out_of_bounds_y = max(seg[:,1]) > (max(y) - 10000)
            plt.plot(seg[:,0], seg[:,1], '.-')
            area = 0.5*np.abs(np.dot(seg[:,0],np.roll(seg[:,1],1))-np.dot(seg[:,1],np.roll(seg[:,0],1))) 
            p = path.Path(seg)
            mask = p.contains_points(X)
            num_points = X[mask].shape[0]
            if num_points >= min_points:
                    seg_list.append([j, area, ii])
            if non_single_cells_x or non_single_cells_y or out_of_bounds_x or out_of_bounds_y or len(seg[:,0])<10:
                continue 
            if area > max_area:
                max_area = area
                best_seg = seg  
                point_num = num_points           
    if best_seg is None or point_num < min_points:
        if len(seg_list) == 0:
            return [samplename, 0, None, 0]
        newbest = None
        smallest_area = np.inf
        for item in seg_list:
            if item[1] < smallest_area:
                smallest_area = item[1]
                best_seg = alls[item[0]][item[2]]
    #fitting an elipse to our identified best fit contour
    vertices1 = np.round(fitEllipse(best_seg[:,0],best_seg[:,1]), 0) 
    #generates the .xml file 
    filename = join(dest, 'gates.xml')
    gate_writer(vertices1, None, None, filename)
    g_strat = fk.parse_gating_xml(filename)    
    gs_results = g_strat.gate_sample(sample)
    #gets the indices of selected cells to move into gate 2
    a = gs_results.get_gate_membership('Polygon1')
    b = Y[a]
    #len= length; a contingency
    if len(b) == 0:
        return [samplename, 0, None, 0]
    plt.scatter(X[:, 0], X[:, 1], c=a, s=12.5)
    ax = plt.gca()
    ax.set_xlabel('SSC-A', fontsize=36)
    ax.set_ylabel('FSC-A', fontsize=36)
    ax.set_title("First Gate", fontsize=30)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=36)
    plt.setp(ax.get_yticklabels(), fontsize=36)
    plt.savefig(join(dest, 'first_gate_KDE.pdf'), format='pdf', bbox_inches='tight')            
    plt.cla()
    plt.clf()
    coefficients = np.polyfit(b[:, 0], b[:, 1], 1)
    poly = np.poly1d(coefficients)
    new_x = np.linspace(min(b[:, 0]), max(b[:, 0]), 2)
    new_y = poly(new_x)
    norms = []
    p1 = np.array([new_x[0], new_y[0]])
    p2 = np.array([new_x[1], new_y[1]])
    for point in b:
        d = norm(np.cross(p2-p1, p1-point))/norm(p2-p1)
        norms.append(d)
    std = np.array(norms).std()
    mask2 = []
    for i in range(len(b)):
        if norms[i] > 4*std:
            mask2.append(False)
        else:
            mask2.append(True)
    coefficients2 = np.polyfit(b[:, 0], b[:, 1], 1)
    poly2 = np.poly1d(coefficients)
    new_x2 = np.linspace(min(b[:, 0]), max(b[:, 0]), 2)
    new_y2 = poly(new_x)   
    factor = 4*std
    plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k')
    plt.scatter(b[:, 0], b[:, 1], c=mask2, s=12.5)
    plt.plot(new_x.tolist(), new_y.tolist(), marker = "o", c='red')
    plt.plot(new_x, [new_y[0] - factor, new_y[1] - factor], marker = "o", c='black')
    plt.plot(new_x, [new_y[0] + factor, new_y[1] + factor], marker = "o", c='black')
    ax = plt.gca()
    ax.set_xlabel('FSC-A', fontsize=36)
    ax.set_ylabel('FSC-H', fontsize=36)
    ax.set_title("Second Gate", fontsize=36)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=36)
    plt.setp(ax.get_yticklabels(), fontsize=36)
    plt.savefig(join(dest, 'secondgate.pdf'), format='pdf', bbox_inches='tight')
    plt.cla()
    plt.clf()
    vertices2 = np.array([[new_x[0], new_y[0] - factor],
            [new_x[1], new_y[1] - factor],
            [new_x[1], new_y[1] + factor],
            [new_x[0], new_y[0] + factor]])
    if np.isnan(vertices1).any() or np.isnan(vertices2).any():
        return [samplename, 0, None, 0]
    gate_writer(vertices1, vertices2, None, filename)
    g_strat = fk.parse_gating_xml(filename)    
    gs_results = g_strat.gate_sample(sample)
    a = gs_results.get_gate_membership('And1')
    c = Y[a]
    fluorescent_chan_names = list(sample.channels['pns'].unique())
    if len(sample.channels) == 6 and 'RFP670' in fluorescent_chan_names:
        return z(samplename, Z, a, vertices1, vertices2, dest, sample, 'RFP670', 'mEmerald', '5-A', '1-A')
    elif len(sample.channels) == 6 and 'mApple' in fluorescent_chan_names:
        return z(samplename, Z, a, vertices1, vertices2, dest, sample, 'mApple', 'mEmerald', 'PE-A', 'FITC-A')
    elif len(sample.channels) == 5:
        return z(samplename, Z, a, vertices1, vertices2, dest, sample, 'mCherry', 'mEmerald', 'PE-Texas Red-A', 'FITC-A')
    elif len(sample.channels) == 7:
        first = z(samplename, Z_ea, a, vertices1, vertices2, dest, sample, 'mApple', 'mEmerald', 'PE-A', 'FITC-A')
        second = z(samplename, Z_er, a, vertices1, vertices2, dest, sample, 'RFP670', 'Emerald', 'APC-A', 'FITC-A')
        third = z(samplename, Z_ar, a, vertices1, vertices2, dest, sample, 'RFP670', 'mApple', 'APC-A', 'PE-A')
        return [first, second, third]


def summarize(outputs, fcs_folder, project_name, skip_renaming):
    dataf = pd.DataFrame(outputs)
    dataf = dataf.rename({0: "File", 1: "FETCH score", 2 : "r_g", 3:"n_tot"}, axis='columns')
    dataf['Dubious?'] = [False for i in range(dataf.shape[0])]
    dataf['FETCH score'] = dataf['FETCH score'].astype(float)
    dataf['r_g'] = dataf['r_g'].astype(float)
    dataf['n_tot'] = dataf['n_tot'].astype(int)
    dataf['Dubious?'] = dataf.apply(lambda row : 'yes' if ((row['r_g'] >= 2) or (row['r_g'] <= 0.5) or (row['n_tot'] < 500) or np.isnan(row['r_g'])) else 'no',
                        axis=1)
    dataf = dataf.drop(['r_g'], axis=1)
    try:
        dataf["File"] = dataf["File"].apply(lambda x: x.rsplit('_')[2] + '_' + x.rsplit('_')[0] + '_' + x.rsplit('_')[1] + '_' + x.rsplit('_')[3] if x not in skip_renaming else x)
    except IndexError:
        pass
    dataf = dataf.sort_values(by='File')
    dataf['Numbername'] = [i for i in range(dataf.shape[0])]
    sns.set(font_scale=2)

    figure(num=None, figsize=(32, 16), dpi=80, facecolor='w', edgecolor='k')
    colors = [(0, 0, 0) if dataf["Dubious?"].iloc[i] == 'no' else (1, 0, 0) for i in range(dataf["File"].unique().shape[0])]

    sns.set_style('ticks')
    sns.catplot(x='File', y='FETCH score', palette=colors, capsize=.2, kind="point", ci="sd", data=dataf, height=10, aspect=2)
    g = sns.swarmplot(x='File', y='FETCH score', data=dataf, color="purple", size=5)
    ax = plt.gca()
    ax.set_xlabel('FETCH_id')
    ax.set_ylabel('FETCH Score')
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=26)
    ax.set(facecolor = "white")
    ax.set_title(project_name, fontsize=26)
    dataf['FETCH score'] = dataf['FETCH score'].fillna(0)
    ax.set(ylim=(0, max(dataf['FETCH score'])+0.1))
    plt.savefig(join(fcs_folder, project_name + ".pdf"), format='pdf', bbox_inches='tight')
    plt.cla()
    plt.clf()
    plt.close()

    dataf = dataf.set_index("File")
    dataf.to_csv(join(fcs_folder, project_name + ".csv"))
#identify which folder contains our fcs files, etc.
def main(args):
    fcs_folder = args.folder
    project_name = args.project
    skip_renaming = args.skip_renaming
    if skip_renaming == '':
        skip_renaming = []
    plt.cla()
    plt.clf()
    plt.close()
    plt.style.use('default')
    inpts = [[join(fcs_folder, samplename), samplename, join(fcs_folder, samplename.rsplit('.')[0])] if (samplename != '.DS_Store' and not os.path.isdir(join(fcs_folder, samplename.rsplit('.')[0]))) else None for samplename in os.listdir(fcs_folder)]
    inpts = list(filter(None, inpts))
    outstuff = []
    for inp in inpts:
        outstuff.append(FETCH_analysis(inp)) 
    outputs = np.array(outstuff)
    #draws the aggregate plot figure and table comparing FETCH scores
    summarize(outputs, fcs_folder, project_name, skip_renaming)
#this is where the code actually starts; runs 'main', above
if __name__ == '__main__':
    args = parse_arguments()
    main(args)
    

