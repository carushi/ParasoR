#!#-*- coding:utf-8 -*-
import sys
import math

""" Minimum value for error"""
eps = math.pow(10, -6)

""" Calculate difference of stem probability between ParasoR with single core and multiple cores."""
def test_diff(single, multiple):
	f = open(single)
	stem_s = [float(line.split('\t')[3]) for line in f.readlines() if len(line) > 0 and line[0] == "*"]
	f.close()
	f = open(multiple)
	stem_m = [float(line.split('\t')[3]) for line in f.readlines() if len(line) > 0 and line[0] == "*"]
	f.close()
	diff = max([abs(i-j) for i, j in zip(stem_s, stem_m)])
	if eps < diff:
		print("Error: Difference between ParasoR with single core and multiple cores is huge! ", diff)
		sys.exit(1)
	else:
		print("Maximum numerical error is under our threshold", diff)
		return


if __name__=="__main__":
	test_diff("../doc/pre.txt", "../doc/para_mem.txt")
	test_diff("../doc/pre_r.txt", "../doc/para_r.txt")
