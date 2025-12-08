# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:23:53 2024

@author: Linoy
"""
import os

ARM_COORDINATES = [1, 124, 250, 344, 493, 584, 692, 742, 883, 932, 1065, 1125,
                   1236, 1296, 1396, 1441, 1542, 1585, 1681, 1721, 1815, 1868,
                   1951, 1987, 2085, 2103, 2200, 2217, 2308, 2327, 2410, 2447,
                   2501, 2527, 2585, 2604, 2666, 2692, 2725, 2753, 2790, 2802,
                   2837, 2852, 2888, 2949, 3045]
chr_arm_start = [1, 251, 495, 694, 886, 1067, 1239, 1399, 1546, 1688, 1824, 1960, 2094, 2210, 2318, 2421, 2512, 2594, 2673, 2733, 2797, 2846, 2898]
ARM_LIST = ['1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q', '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q', '11p', '11q', '12p', '12q', '13p', '13q', '14p', '14q', 
            '15p', '15q', '16p', '16q', '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21p', '21q', '22p', '22q', 'Xp', 'Xq']
PERCENT_CUTOFF = 75


CUSTOM_COLORS = {
    "Colorectal": "#6268F0",
    "Breast": "#F0A8CC",
    "gastric": "#a6b7e3",
    "Lung": "#FFFFFF",
    "bladder": "#e0c02f",
    "Head and Neck": "#E05A57",
    "sarcoma": "#e4e80c",
    "Glioblastoma": "#BEBEBE",
    "Melanoma": "#808080",
    "HCC": "#82d7b4",
    "Pancreatic": "#B97AF0",
    "Renal": "#E1AF7D",
    "Esophageal": "#66E0DB",
    "Epithelial":"#66cdaa",
    "Prostate":'#8CAAEF',
    "Uterine": "#FFDAB9",
    "Leiomyosarcoma": "#800080",
    "Liver":"#50C878",
    "Neuroendocrine": "#800080",  
    "Ovarian": "#008080",
    "Bladder": "#FFFEA0",
    "Gastric": "#a6b7e3",
    "Salivary": "#CCCCFF",
    "Skin": "#000000"}


