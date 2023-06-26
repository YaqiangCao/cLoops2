#!/usr/bin/env python
#--coding:utf-8 --
"""
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import warnings
warnings.filterwarnings("ignore")

#3rd plotting setting
import matplotlib as mpl
mpl.use("pdf")
import seaborn as sns
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 8.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import pylab
sns.set_style("white")
from matplotlib.colors import ListedColormap

#colors defined from brewer2mpl
colors = [
    (0.8941176470588236, 0.10196078431372549, 0.10980392156862745),
    (0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
    (0.30196078431372547, 0.6862745098039216, 0.2901960784313726),
    (0.596078431372549, 0.3058823529411765, 0.6392156862745098),
    (1.0, 0.4980392156862745, 0.0),
    (0,0,0), #black
    #(1.0, 1.0, 0.2),
    #(0.7, 1.0, 0.2),
    (0.6509803921568628, 0.33725490196078434, 0.1568627450980392),
    (0.9686274509803922, 0.5058823529411764, 0.7490196078431373),
    (0.6, 0.6, 0.6),
    (0.5529411764705883, 0.8274509803921568, 0.7803921568627451),
    (1.0, 1.0, 0.7019607843137254),
    (0.7450980392156863, 0.7294117647058823, 0.8549019607843137),
    (0.984313725490196, 0.5019607843137255, 0.4470588235294118),
    (0.5019607843137255, 0.6941176470588235, 0.8274509803921568),
    (0.9921568627450981, 0.7058823529411765, 0.3843137254901961),
    (0.7019607843137254, 0.8705882352941177, 0.4117647058823529),
    (0.9882352941176471, 0.803921568627451, 0.8980392156862745),

]

cmap1 = sns.diverging_palette(250, 15, s=75, l=40, n=9).as_hex()
cmap1[int(len(cmap1) / 2)] = "#FFFFFF"
cmap1 = ListedColormap(cmap1)
cmap2 = sns.color_palette("RdBu_r", 11).as_hex()
cmap2[int(len(cmap2) / 2)] = "#FFFFFF"
cmap2 = ListedColormap(cmap2)
cmap3 = sns.color_palette("coolwarm", 11).as_hex()
cmap3[int(len(cmap3) / 2)] = "#FFFFFF"
cmap3 = ListedColormap(cmap3)
cmaps = {
    "red": sns.cubehelix_palette(light=1, as_cmap=True),
    "div": cmap1,
    "summer": cmap2,
    "cool": cmap3,
}
