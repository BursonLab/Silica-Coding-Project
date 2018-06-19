#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 12:56:06 2017

@author: DanWall
"""


def xyx_to_list(file):
    text = []
    file_lines = file.readlines()
    file_lines = [x.strip() for x in file_lines]
    file_lines = file_lines[2:]
    for line in file_lines:
        row = []
        rowNum = int(line[0])
        row.append(rowNum)
        x = float(line[5:14])
        y = float(line[18:27])
        z = float(line[31:40])
        row.append(x)
        row.append(y)
        row.append(z)
        text.append(row)
    return text


def main():
    text = xyx_to_list(open('Sample xyz template.txt', encoding = 'utf-8'))
    for row in text:
        print(row)


main()
