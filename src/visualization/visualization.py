import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import json
import subprocess

with open("config.json") as f:
    cfg = json.load(f)

draw_survival_plot = cfg['draw_survival_plot']
draw_focal_plots = cfg['draw_focal_score_plots']
draw_focal_proportion_plots = cfg['draw_focal_score_probability_plots']

if draw_survival_plot:
    subprocess.call("/src/visualization/survival_plots.R", shell=True)

lgg_p = 'data/external/LGG.focal_score_by_genes.txt'
gbm_p = 'data/external/GBM.focal_score_by_genes.txt'
lgg_df = pd.read_csv(lgg_p, sep = '\t')
gbm_df = pd.read_csv(gbm_p, sep = '\t')


if draw_focal_plots:
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist = [(1, 0, 0, 1), (0, 0, 0, 0), (0, 0, 1, 1)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)

    plt.imshow(lgg_df.iloc[:, 3:].T, cmap=cmap, aspect='auto')
    plt.title('LGG focal score by gene')
    plt.xlabel('Gene Indices')
    plt.ylabel('Case Indices')
    plt.savefig('reports/figures/LGG focal map.png', dpi=300)

    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist = [(1, 0, 0, 1), (0, 0, 0, 0), (0, 0, 1, 1)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)

    plt.imshow(gbm_df.iloc[:, 3:].T, cmap=cmap, aspect='auto')
    plt.title('GBM Focal Score by Gene')
    plt.xlabel('Gene Indices')
    plt.ylabel('Case Indices')
    plt.savefig('reports/figures/GBM focal map.png', dpi=300)
    
if draw_focal_proportion_plots:
    b, indices = set(), []
    for i, e in enumerate(gbm_df['Cytoband']):
        e = e.replace('p', 'q').split('q')[0]
        if e not in b:
            b.add(e)
            indices.append(i)

    scale = .022
    gbm_x = gbm_df.iloc[:, 3:].T.columns
    gbm_dup = gbm_df.iloc[:, 3:].T.apply(lambda c: len(c[c>0]))
    gbm_del = gbm_df.iloc[:, 3:].T.apply(lambda c: len(c[c<0])) * -1
    gbm_dup = gbm_dup / len(gbm_df)
    gbm_del = gbm_del / len(gbm_df)
    
    plt.scatter(gbm_x, gbm_dup, s=.2, color='blue')
    plt.scatter(gbm_x, gbm_del, s=.2, color='red')
    gbm_axes = plt.gca()
    gbm_axes.set_ylim([scale*-1, scale])
    plt.title('GBM Average Probabilities of Duplications and Deletions by Gene Indices')
    plt.xlabel('Gene Indices')
    plt.ylabel("Focal Score Proportions")
    for i in indices:
        plt.axvline(x=i, color='black', linewidth=.5)
    plt.savefig('reports/figures/GBM probabilities zoomed and delineated', dpi=300)
    
    
    b, indices = set(), []
    for i, e in enumerate(lgg_df['Cytoband']):
        e = e.replace('p', 'q').split('q')[0]
        if e not in b:
            b.add(e)
            indices.append(i)
            
    lgg_x = lgg_df.iloc[:, 3:].T.columns
    lgg_dup = lgg_df.iloc[:, 3:].T.apply(lambda c: len(c[c>0]))
    lgg_del = lgg_df.iloc[:, 3:].T.apply(lambda c: len(c[c<0])) * -1
    lgg_dup = lgg_dup / len(lgg_df)
    lgg_del = lgg_del / len(lgg_df)
    
    lgg_axes = plt.gca()
    lgg_axes.set_ylim([scale*-1, scale])

    plt.scatter(lgg_x, lgg_dup, s=.2, color='blue')
    plt.scatter(lgg_x, lgg_del, s=.2, color='red')
    plt.title('LGG Average Probabilities of Duplications and Deletions by Gene Indices')
    plt.xlabel('Gene Indices')
    plt.ylabel("Focal Score Proportions")
    for i in indices:
        plt.axvline(x=i, color='black', linewidth=.5)
    plt.savefig('reports/figures/LGG probabilities zoomed and delineated', dpi=300)