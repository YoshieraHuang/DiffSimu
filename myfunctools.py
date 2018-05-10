# -*- coding: utf-8 -*-
# @Date    : 2017-12-27 20:45:56
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np


def where(arr, value):
	idx = np.argwhere(arr == value)
	if idx.shape == (0,1):
		return 0
	return idx[0,0] + 1

def equal(v1, v2, EPS):
	if np.fabs(v1 - v2) > EPS:
		return False
		
	return True

if __name__ == '__main__':
	pass