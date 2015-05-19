#!/usr/bin/env python


import argparse
from pyworkflow.em.plotter import plotFile
from pyworkflow.gui.plotter import Plotter

def main():
    parser = argparse.ArgumentParser(prog='Scipion Plot')
    parser.add_argument('--file', help='File to visualize', required=True)
    parser.add_argument('--block', help='Block to visualize')
    parser.add_argument('--type', help='Plot type')
    parser.add_argument('--columns', help='Columns to plot')
    parser.add_argument('--xcolumn', help='X Column to plot')
    parser.add_argument('--orderColumn', help='Column to order')
    parser.add_argument('--orderDir', help='Order direction(ASC, DESC)')
    parser.add_argument('--bins', help='If plot type is histogram, number of bins')
    parser.add_argument('--colors', help='Colors to plot columns')
    parser.add_argument('--styles', help='Styles to plot columns')
    parser.add_argument('--markers', help='Markers to plot columns')
    parser.add_argument('--title', help='Plot title', default='')
    parser.add_argument('--ytitle', help='Y axis title', default='')
    parser.add_argument('--xtitle', help='X axis title', default='')


    
    args = parser.parse_args()
    #print args
    plotfile = args.file
    block = args.block if args.block else '' 
    type = args.type
    columns = args.columns
    xcolumn = args.xcolumn
    orderColumn = args.orderColumn
    orderDir = args.orderDir
    
    bins = args.bins
    colors = args.colors
    styles = args.styles
    markers = args.markers
    title = args.title
    xtitle = args.xtitle
    ytitle = args.ytitle
    
    Plotter.setBackend('TkAgg')
    plotFile(plotfile, block, type,
               columns, colors, styles, markers,
               xcolumn, ytitle, xtitle, title, bins, orderColumn, orderDir).show(block=True)
#     else: 
#         plotMetaData(plotfile, *args)
    
    
    
if __name__ == '__main__':
    main()
    
    

