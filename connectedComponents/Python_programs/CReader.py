#!/usr/bin/env python
# -*- coding: utf-8 -*- import sys

import csv, json

def readCsvToArray(csv_file):
	return readCsvToArrayDelimiter(csv_file, ",")
def readCsvHeader(csv_file):
	return readCsvToArrayDelimiter(csv_file, ",")
def readCsvToDic(csv_file):
	return csvArrayToDic(readCsvToArray(csv_file))
def readCsvToDicDelimiter(csv_file, d):
	return csvArrayToDic(readCsvToArrayDelimiter(csv_file, d))
def readCsvToDicHeader(csv_file, header):
	return csvArrayToDic(csvArrayChangeHeader(readCsvToArray(csv_file), header))
def readCsvToDicDelimiterHeader(csv_file, d, header):
	return csvArrayToDic(csvArrayChangeHeader(readCsvToArrayDelimiter(csv_file, d), header))
def csvArrayStringToIntByKey(rows, keys):
	return csvArrayStringToIntByIndex(rows, getIndexFromKey(rows[0], keys))
def csvArrayStringToFloatByKey(rows, keys):
	return csvArrayStringToFloatByIndex(rows, getIndexFromKey(rows[0], keys))
def csvArrayToDic(rows):
	return csvArrayToDicHeader(rows, rows[0])
def readCsvToArrayDelimiter(csv_file, d):
    data = []
    with open(csv_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=d)
        for row in readCSV:
            data.append(row)
    return data
def readCsvHeader(csv_file):
    data = []
    with open(csv_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=",")
        for row in readCSV:
        	return row
    return data

def csvArrayChangeHeader(rows, header):
	rows[0] = header
	return rows
def getIndexFromKey(header, keys):
	n = len(header)
	keyMap = {}
	for i in range(n):
		keyMap[header[i]] = i
	indexes = []
	for k in keys:
		indexes.append(keyMap[k])
	return indexes
def csvArrayStringToIntByIndex(rows, indexes):
	n = len(rows)
	for i in range(1, n):
		for j in indexes:
			rows[i][j] = int(rows[i][j])
	return rows
def csvArrayStringToFloatByIndex(rows, indexes):
	n = len(rows)
	for i in range(1, n):
		for j in indexes:
			rows[i][j] = float(rows[i][j])
	return rows
def csvArrayToDicHeader(rows, header):
	data = []
	n = len(rows)
	h = len(header)
	for i in range(1,n):
		v = {}
		for j in range(h):
			v[header[j]] = rows[i][j]
		data.append(v)
	return data
def dicStrToFloat(data, keys):
	if len(rows[0]) != len(rows[1]): print ("###csvArrayToDicByKeys### heder is wrong.")
	header = rows[0]
	n = len(header)
	tags = {}
	for i in range(n): tags[header[i]] = False
	for k in floatTags: 
		if k in tags: tags[k] = True
	data = []
	for i in range(1, len(rows)):
		v = {}
		for j in range(n):
			k = header[j]
			u = rows[i][j]
			if tags[k] == True: 
				if type(u) == int: v[k] = int(u) 
				else: v[k] = float(u) 
			else: v[k] = u
		data.append(v)
	return data
def readJson(json_file):
	json_data = open(json_file)
	values = json.load(json_data)
	json_data.close()
	print ("DONE load ", json_file)
	print 
	return values
def saveJson(json_file, data):
	with open(json_file, 'w') as outfile:
		json.dump(data, outfile)	
	print ("DONE save ", json_file)
	print 
def saveCsv(csv_file, data):
	with open(csv_file, 'w', newline='') as outfile:
	    writer = csv.writer(outfile)
	    for row in data:
	    	writer.writerow(row)