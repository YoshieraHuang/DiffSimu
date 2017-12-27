# -*- coding: utf-8 -*-
# @Date    : 2017-12-27 20:45:56
# @Author  : J. W. Huang (huangjasper@126.com)

def tolist(gen):
	'''
		list all elements in generator
	'''
	if not type(gen) is generator:
		raise TypeError('Argument must be a generator')

	results = list(gen);
	if len(results) == 1:
		return results[0]
	else:
		return results


if __name__ == '__main__':
	pass