import numpy as np

model='xpt_nnlo_FV'

dr = 0.5
p_range = np.arange(.1,10+dr,dr)
a_range = np.arange(.1,10+dr,dr)

nnlo_width = dict()

nnlo_width['xpt_nnlo_FV_alphaS'] = dict()
nnlo_width['xpt_nnlo_FV_alphaS']['PK'] = dict()
nnlo_width['xpt_nnlo_FV_alphaS']['PK']['p_4']   = 3.5
nnlo_width['xpt_nnlo_FV_alphaS']['PK']['k_4']   = 3.5
nnlo_width['xpt_nnlo_FV_alphaS']['PK']['s_4']   = 2.7
nnlo_width['xpt_nnlo_FV_alphaS']['PK']['saS_4'] = 2.7

nnlo_width['xpt_nnlo_FV_alphaS_logSq'] = dict()
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK'] = dict()
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK']['p_4']   = 1.9
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK']['k_4']   = 1.9
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK']['s_4']   = 2.9
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK']['saS_4'] = 2.9

nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq'] = dict()
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK'] = dict()
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK']['p_4']   = 9.8
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK']['k_4']   = 9.8
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK']['s_4']   = 3.3
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK']['saS_4'] = 3.3

nnlo_width['ma_nnlo_FV_alphaS_logSq'] = dict()
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK'] = dict()
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK']['p_4']   = 0.9
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK']['k_4']   = 0.9
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK']['s_4']   = 0.8
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK']['saS_4'] = 0.8

nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq'] = dict()
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK'] = dict()
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK']['p_4']   = 0.9
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK']['k_4']   = 0.9
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK']['s_4']   = 0.1
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK']['saS_4'] = 0.1
