#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 12:56:06 2017

@author: DanWall
"""
import Si_Ring_Classes

def list_of_centers_to_objects(center_xyz_list):
    center_obj_list = []
    for c in center_xyz_list:
        center = Si_Ring_Classes.ring_center(c[0], c[1], c[2], c[3])
        center_obj_list.append(center)
        

def xyz_to_list(file):
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
    text = xyz_to_list(open('Sample xyz template.txt', encoding = 'utf-8'))
    for row in text:
        print(row)


main()
