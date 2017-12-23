# -*- coding: utf-8 -*-
# @Date    : 2017-12-23 13:48:19
# @Author  : J. W. Huang (huangjasper@126.com)


dict_material = {
	'Cu' : ('fcc', (3.615, 'Cu'))
}

dict_element = {
	'Cu' : ('Cu', 64, 231),
	'Fe' : ('Fe', 54, 215)
}

dict_structure = {
	'fcc' : (('$a', '$a', '$a', 90, 90, 90), (4, ('$elem', (0,0,0)), ('$elem', (0.5,0.5,0)), ('$elem', (0.5,0,0.5)), ('$elem', (0,0.5,0.5)))),
	'bcc' : (('$a', '$a', '$a', 90, 90, 90), (2, ('$elem', (0,0,0)), ('$elem', (0.5,0.5,0.5)))),
	'sc'  : (('$a', '$a', '$a', 90, 90, 90), (1, ('$elem', (0,0,0)))),
}