# -*- coding: utf-8 -*-
# @Date    : 2017-12-27 20:45:56
# @Author  : J. W. Huang (huangjasper@126.com)

import types
import numpy as np

def tolist(gen):

	if not isinstance(gen, types.GeneratorType):
		raise TypeError('Argument must be a generator')

	l = list(gen)
	return l[0] if len(l) == 1 else l

if __name__ == '__main__':
	pass