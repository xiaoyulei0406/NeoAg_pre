import os, sys, subprocess
import argparse
import glob

IEDB='/home/cyu/tools/IEDB/mhc_ii/mhc_II_binding.py '
methods=['IEDB_recommended','sturniolo','consensus3','nn_align','smm_align','comblib','netmhciipan_el','netmhciipan_ba']

allele_netmhcII='HLA-DPA1*01:03/DPB1*02:01,HLA-DPA1*01:03/DPB1*03:01,HLA-DPA1*01:03/DPB1*04:01,HLA-DPA1*01:03/DPB1*04:02,HLA-DPA1*01:03/DPB1*06:01,HLA-DPA1*01:03/DPB1*11:01,HLA-DPA1*01:03/DPB1*17:01,HLA-DPA1*01:03/DPB1*20:01,HLA-DPA1*01:03/DPB1*23:01,HLA-DQA1*01:02/DQB1*05:01,HLA-DQA1*01:02/DQB1*06:02,HLA-DQA1*01:02/DQB1*06:04,HLA-DQA1*01:03/DQB1*05:01,HLA-DQA1*01:03/DQB1*06:03,HLA-DQA1*02:01/DQB1*02:01,HLA-DQA1*02:01/DQB1*02:02,HLA-DQA1*04:01/DQB1*03:01,HLA-DQA1*05:01/DQB1*02:01,HLA-DQA1*05:01/DQB1*03:01,HLA-DQA1*05:05/DQB1*03:01,DRB1*01:01,DRB1*01:02,DRB1*01:03,DRB1*01:04,DRB1*01:05,DRB1*01:06,DRB1*01:07,DRB1*01:08,DRB1*01:09,DRB1*01:10,DRB1*01:11,DRB1*01:12,DRB1*01:13,DRB1*01:14,DRB1*01:15,DRB1*01:16,DRB1*01:17,DRB1*01:18,DRB1*01:19,DRB1*01:20,DRB1*01:21,DRB1*01:22,DRB1*01:23,DRB1*01:24,DRB1*01:25,DRB1*01:26,DRB1*01:27,DRB1*01:28,DRB1*01:29,DRB1*01:30,DRB1*01:31,DRB1*01:32,DRB1*03:01,DRB1*03:02,DRB1*03:03,DRB1*03:04,DRB1*03:05,DRB1*03:06,DRB1*03:07,DRB1*03:08,DRB1*03:10,DRB1*03:11,DRB1*03:13,DRB1*03:14,DRB1*03:15,DRB1*03:17,DRB1*03:18,DRB1*03:19,DRB1*03:20,DRB1*03:21,DRB1*03:22,DRB1*03:23,DRB1*03:24,DRB1*03:25,DRB1*03:26,DRB1*03:27,DRB1*03:28,DRB1*03:29,DRB1*03:30,DRB1*03:31,DRB1*03:32,DRB1*03:33,DRB1*03:34,DRB1*03:35,DRB1*03:36,DRB1*03:37,DRB1*03:38,DRB1*03:39,DRB1*03:40,DRB1*03:41,DRB1*03:42,DRB1*03:43,DRB1*03:44,DRB1*03:45,DRB1*03:46,DRB1*03:47,DRB1*03:48,DRB1*03:49,DRB1*03:50,DRB1*03:51,DRB1*03:52,DRB1*03:53,DRB1*03:54,DRB1*03:55,DRB1*04:01,DRB1*04:02,DRB1*04:03,DRB1*04:04,DRB1*04:05,DRB1*04:06,DRB1*04:07,DRB1*04:08,DRB1*04:09,DRB1*04:10,DRB1*04:11,DRB1*04:12,DRB1*04:13,DRB1*04:14,DRB1*04:15,DRB1*04:16,DRB1*04:17,DRB1*04:18,DRB1*04:19,DRB1*04:21,DRB1*04:22,DRB1*04:23,DRB1*04:24,DRB1*04:26,DRB1*04:27,DRB1*04:28,DRB1*04:29,DRB1*04:30,DRB1*04:31,DRB1*04:33,DRB1*04:34,DRB1*04:35,DRB1*04:36,DRB1*04:37,DRB1*04:38,DRB1*04:39,DRB1*04:40,DRB1*04:41,DRB1*04:42,DRB1*04:43,DRB1*04:44,DRB1*04:45,DRB1*04:46,DRB1*04:47,DRB1*04:48,DRB1*04:49,DRB1*04:50,DRB1*04:51,DRB1*04:52,DRB1*04:53,DRB1*04:54,DRB1*04:55,DRB1*04:56,DRB1*04:57,DRB1*04:58,DRB1*04:59,DRB1*04:60,DRB1*04:61,DRB1*04:62,DRB1*04:63,DRB1*04:64,DRB1*04:65,DRB1*04:66,DRB1*04:67,DRB1*04:68,DRB1*04:69,DRB1*04:70,DRB1*04:71,DRB1*04:72,DRB1*04:73,DRB1*04:74,DRB1*04:75,DRB1*04:76,DRB1*04:77,DRB1*04:78,DRB1*04:79,DRB1*04:80,DRB1*04:82,DRB1*04:83,DRB1*04:84,DRB1*04:85,DRB1*04:86,DRB1*04:87,DRB1*04:88,DRB1*04:89,DRB1*04:91,DRB1*07:01,DRB1*07:03,DRB1*07:04,DRB1*07:05,DRB1*07:06,DRB1*07:07,DRB1*07:08,DRB1*07:09,DRB1*07:11,DRB1*07:12,DRB1*07:13,DRB1*07:14,DRB1*07:15,DRB1*07:16,DRB1*07:17,DRB1*07:19,DRB1*08:01,DRB1*08:02,DRB1*08:03,DRB1*08:04,DRB1*08:05,DRB1*08:06,DRB1*08:07,DRB1*08:08,DRB1*08:09,DRB1*08:10,DRB1*08:11,DRB1*08:12,DRB1*08:13,DRB1*08:14,DRB1*08:15,DRB1*08:16,DRB1*08:18,DRB1*08:19,DRB1*08:20,DRB1*08:21,DRB1*08:22,DRB1*08:23,DRB1*08:24,DRB1*08:25,DRB1*08:26,DRB1*08:27,DRB1*08:28,DRB1*08:29,DRB1*08:30,DRB1*08:31,DRB1*08:32,DRB1*08:33,DRB1*08:34,DRB1*08:35,DRB1*08:36,DRB1*08:37,DRB1*08:38,DRB1*08:39,DRB1*08:40,DRB1*09:01,DRB1*09:02,DRB1*09:03,DRB1*09:04,DRB1*09:05,DRB1*09:06,DRB1*09:07,DRB1*09:08,DRB1*09:09,DRB1*10:01,DRB1*10:02,DRB1*10:03,DRB1*11:01,DRB1*11:02,DRB1*11:03,DRB1*11:04,DRB1*11:05,DRB1*11:06,DRB1*11:07,DRB1*11:08,DRB1*11:09,DRB1*11:10,DRB1*11:11,DRB1*11:12,DRB1*11:13,DRB1*11:14,DRB1*11:15,DRB1*11:16,DRB1*11:17,DRB1*11:18,DRB1*11:19,DRB1*11:20,DRB1*11:21,DRB1*11:24,DRB1*11:25,DRB1*11:27,DRB1*11:28,DRB1*11:29,DRB1*11:30,DRB1*11:31,DRB1*11:32,DRB1*11:33,DRB1*11:34,DRB1*11:35,DRB1*11:36,DRB1*11:37,DRB1*11:38,DRB1*11:39,DRB1*11:41,DRB1*11:42,DRB1*11:43,DRB1*11:44,DRB1*11:45,DRB1*11:46,DRB1*11:47,DRB1*11:48,DRB1*11:49,DRB1*11:50,DRB1*11:51,DRB1*11:52,DRB1*11:53,DRB1*11:54,DRB1*11:55,DRB1*11:56,DRB1*11:57,DRB1*11:58,DRB1*11:59,DRB1*11:60,DRB1*11:61,DRB1*11:62,DRB1*11:63,DRB1*11:64,DRB1*11:65,DRB1*11:66,DRB1*11:67,DRB1*11:68,DRB1*11:69,DRB1*11:70,DRB1*11:72,DRB1*11:73,DRB1*11:74,DRB1*11:75,DRB1*11:76,DRB1*11:77,DRB1*11:78,DRB1*11:79,DRB1*11:80,DRB1*11:81,DRB1*11:82,DRB1*11:83,DRB1*11:84,DRB1*11:85,DRB1*11:86,DRB1*11:87,DRB1*11:88,DRB1*11:89,DRB1*11:90,DRB1*11:91,DRB1*11:92,DRB1*11:93,DRB1*11:94,DRB1*11:95,DRB1*11:96,DRB1*12:01,DRB1*12:02,DRB1*12:03,DRB1*12:04,DRB1*12:05,DRB1*12:06,DRB1*12:07,DRB1*12:08,DRB1*12:09,DRB1*12:10,DRB1*12:11,DRB1*12:12,DRB1*12:13,DRB1*12:14,DRB1*12:15,DRB1*12:16,DRB1*12:17,DRB1*12:18,DRB1*12:19,DRB1*12:20,DRB1*12:21,DRB1*12:22,DRB1*12:23,DRB1*13:01,DRB1*13:02,DRB1*13:03,DRB1*13:04,DRB1*13:05,DRB1*13:06,DRB1*13:07,DRB1*13:08,DRB1*13:09,DRB1*13:10,DRB1*13:100,DRB1*13:101,DRB1*13:11,DRB1*13:12,DRB1*13:13,DRB1*13:14,DRB1*13:15,DRB1*13:16,DRB1*13:17,DRB1*13:18,DRB1*13:19,DRB1*13:20,DRB1*13:21,DRB1*13:22,DRB1*13:23,DRB1*13:24,DRB1*13:26,DRB1*13:27,DRB1*13:29,DRB1*13:30,DRB1*13:31,DRB1*13:32,DRB1*13:33,DRB1*13:34,DRB1*13:35,DRB1*13:36,DRB1*13:37,DRB1*13:38,DRB1*13:39,DRB1*13:41,DRB1*13:42,DRB1*13:43,DRB1*13:44,DRB1*13:46,DRB1*13:47,DRB1*13:48,DRB1*13:49,DRB1*13:50,DRB1*13:51,DRB1*13:52,DRB1*13:53,DRB1*13:54,DRB1*13:55,DRB1*13:56,DRB1*13:57,DRB1*13:58,DRB1*13:59,DRB1*13:60,DRB1*13:61,DRB1*13:62,DRB1*13:63,DRB1*13:64,DRB1*13:65,DRB1*13:66,DRB1*13:67,DRB1*13:68,DRB1*13:69,DRB1*13:70,DRB1*13:71,DRB1*13:72,DRB1*13:73,DRB1*13:74,DRB1*13:75,DRB1*13:76,DRB1*13:77,DRB1*13:78,DRB1*13:79,DRB1*13:80,DRB1*13:81,DRB1*13:82,DRB1*13:83,DRB1*13:84,DRB1*13:85,DRB1*13:86,DRB1*13:87,DRB1*13:88,DRB1*13:89,DRB1*13:90,DRB1*13:91,DRB1*13:92,DRB1*13:93,DRB1*13:94,DRB1*13:95,DRB1*13:96,DRB1*13:97,DRB1*13:98,DRB1*13:99,DRB1*14:01,DRB1*14:02,DRB1*14:03,DRB1*14:04,DRB1*14:05,DRB1*14:06,DRB1*14:07,DRB1*14:08,DRB1*14:09,DRB1*14:10,DRB1*14:11,DRB1*14:12,DRB1*14:13,DRB1*14:14,DRB1*14:15,DRB1*14:16,DRB1*14:17,DRB1*14:18,DRB1*14:19,DRB1*14:20,DRB1*14:21,DRB1*14:22,DRB1*14:23,DRB1*14:24,DRB1*14:25,DRB1*14:26,DRB1*14:27,DRB1*14:28,DRB1*14:29,DRB1*14:30,DRB1*14:31,DRB1*14:32,DRB1*14:33,DRB1*14:34,DRB1*14:35,DRB1*14:36,DRB1*14:37,DRB1*14:38,DRB1*14:39,DRB1*14:40,DRB1*14:41,DRB1*14:42,DRB1*14:43,DRB1*14:44,DRB1*14:45,DRB1*14:46,DRB1*14:47,DRB1*14:48,DRB1*14:49,DRB1*14:50,DRB1*14:51,DRB1*14:52,DRB1*14:53,DRB1*14:54,DRB1*14:55,DRB1*14:56,DRB1*14:57,DRB1*14:58,DRB1*14:59,DRB1*14:60,DRB1*14:61,DRB1*14:62,DRB1*14:63,DRB1*14:64,DRB1*14:65,DRB1*14:67,DRB1*14:68,DRB1*14:69,DRB1*14:70,DRB1*14:71,DRB1*14:72,DRB1*14:73,DRB1*14:74,DRB1*14:75,DRB1*14:76,DRB1*14:77,DRB1*14:78,DRB1*14:79,DRB1*14:80,DRB1*14:81,DRB1*14:82,DRB1*14:83,DRB1*14:84,DRB1*14:85,DRB1*14:86,DRB1*14:87,DRB1*14:88,DRB1*14:89,DRB1*14:90,DRB1*14:91,DRB1*14:93,DRB1*14:94,DRB1*14:95,DRB1*14:96,DRB1*14:97,DRB1*14:98,DRB1*14:99,DRB1*15:01,DRB1*15:02,DRB1*15:03,DRB1*15:04,DRB1*15:05,DRB1*15:06,DRB1*15:07,DRB1*15:08,DRB1*15:09,DRB1*15:10,DRB1*15:11,DRB1*15:12,DRB1*15:13,DRB1*15:14,DRB1*15:15,DRB1*15:16,DRB1*15:18,DRB1*15:19,DRB1*15:20,DRB1*15:21,DRB1*15:22,DRB1*15:23,DRB1*15:24,DRB1*15:25,DRB1*15:26,DRB1*15:27,DRB1*15:28,DRB1*15:29,DRB1*15:30,DRB1*15:31,DRB1*15:32,DRB1*15:33,DRB1*15:34,DRB1*15:35,DRB1*15:36,DRB1*15:37,DRB1*15:38,DRB1*15:39,DRB1*15:40,DRB1*15:41,DRB1*15:42,DRB1*15:43,DRB1*15:44,DRB1*15:45,DRB1*15:46,DRB1*15:47,DRB1*15:48,DRB1*15:49,DRB1*16:01,DRB1*16:02,DRB1*16:03,DRB1*16:04,DRB1*16:05,DRB1*16:07,DRB1*16:08,DRB1*16:09,DRB1*16:10,DRB1*16:11,DRB1*16:12,DRB1*16:14,DRB1*16:15,DRB1*16:16,DRB3*01:01,DRB3*01:04,DRB3*01:05,DRB3*01:08,DRB3*01:09,DRB3*01:11,DRB3*01:12,DRB3*01:13,DRB3*01:14,DRB3*02:01,DRB3*02:02,DRB3*02:04,DRB3*02:05,DRB3*02:09,DRB3*02:10,DRB3*02:11,DRB3*02:12,DRB3*02:13,DRB3*02:14,DRB3*02:15,DRB3*02:16,DRB3*02:17,DRB3*02:18,DRB3*02:19,DRB3*02:20,DRB3*02:21,DRB3*02:22,DRB3*02:23,DRB3*02:24,DRB3*02:25,DRB3*03:01,DRB3*03:03,DRB4*01:01,DRB4*01:03,DRB4*01:04,DRB4*01:06,DRB4*01:07,DRB4*01:08,DRB5*01:01,DRB5*01:02,DRB5*01:03,DRB5*01:04,DRB5*01:05,DRB5*01:06,DRB5*01:08N,DRB5*01:11,DRB5*01:12,DRB5*01:13,DRB5*01:14,DRB5*02:02,DRB5*02:03,DRB5*02:04,DRB5*02:05'
allele_sturniolo='HLA-DRB1*01:01,HLA-DRB1*04:21,HLA-DRB1*11:07,HLA-DRB1*13:28,HLA-DRB1*01:02,HLA-DRB1*04:23,HLA-DRB1*11:14,HLA-DRB1*15:01,HLA-DRB1*03:01,HLA-DRB1*04:26,HLA-DRB1*11:20,HLA-DRB1*15:02,HLA-DRB1*03:05,HLA-DRB1*07:01,HLA-DRB1*11:21,HLA-DRB1*15:06,HLA-DRB1*03:06,HLA-DRB1*07:03,HLA-DRB1*11:28,HLA-DRB5*01:01,HLA-DRB1*03:07,HLA-DRB1*08:01,HLA-DRB1*13:01,HLA-DRB5*01:05,HLA-DRB1*03:08,HLA-DRB1*08:02,HLA-DRB1*13:02,HLA-DRB1*03:09,HLA-DRB1*08:04,HLA-DRB1*13:04,HLA-DRB1*03:11,HLA-DRB1*08:06,HLA-DRB1*13:05,HLA-DRB1*04:01,HLA-DRB1*08:13,HLA-DRB1*13:07,HLA-DRB1*04:02,HLA-DRB1*08:17,HLA-DRB1*13:11,HLA-DRB1*04:04,HLA-DRB1*11:01,HLA-DRB1*13:21,HLA-DRB1*04:05,HLA-DRB1*11:02,HLA-DRB1*13:22,HLA-DRB1*04:08,HLA-DRB1*11:04,HLA-DRB1*13:23,HLA-DRB1*04:10,HLA-DRB1*11:06,HLA-DRB1*13:27'
allele_nn_align='HLA-DPA1*01/DPB1*04:01,HLA-DRB1*01:01,HLA-DRB1*15:01,HLA-DPA1*01:03/DPB1*02:01,HLA-DRB1*03:01,HLA-DRB3*01:01,HLA-DPA1*02:01/DPB1*01:01,HLA-DRB1*04:01,HLA-DRB4*01:01,HLA-DPA1*02:01/DPB1*05:01,HLA-DRB1*04:04,HLA-DRB5*01:01,HLA-DPA1*03:01/DPB1*04:02,HLA-DRB1*04:05,HLA-DQA1*01:01/DQB1*05:01,HLA-DRB1*07:01,HLA-DQA1*01:02/DQB1*06:02,HLA-DRB1*08:02,HLA-DQA1*03:01/DQB1*03:02,HLA-DRB1*09:01,HLA-DQA1*04:01/DQB1*04:02,HLA-DRB1*11:01,HLA-DQA1*05:01/DQB1*02:01,HLA-DRB1*12:01,HLA-DQA1*05:01/DQB1*03:01,HLA-DRB1*13:02'
allele_smm_align='HLA-DPA1*01/DPB1*04:01,HLA-DRB1*01:01,HLA-DRB1*15:01,HLA-DPA1*01:03/DPB1*02:01,HLA-DRB1*03:01,HLA-DRB3*01:01,HLA-DPA1*02:01/DPB1*01:01,HLA-DRB1*04:01,HLA-DRB4*01:01,HLA-DPA1*02:01/DPB1*05:01,HLA-DRB1*04:04,HLA-DRB5*01:01,HLA-DPA1*03:01/DPB1*04:02,HLA-DRB1*04:05,HLA-DQA1*01:01/DQB1*05:01,HLA-DRB1*07:01,HLA-DQA1*01:02/DQB1*06:02,HLA-DRB1*08:02,HLA-DQA1*03:01/DQB1*03:02,HLA-DRB1*09:01,HLA-DQA1*04:01/DQB1*04:02,HLA-DRB1*11:01,HLA-DQA1*05:01/DQB1*02:01,HLA-DRB1*12:01,HLA-DQA1*05:01/DQB1*03:01,HLA-DRB1*13:02'
allele_comblib='HLA-DPA1*01/DPB1*04:01,HLA-DRB1*01:01,HLA-DPA1*01:03/DPB1*02:01,HLA-DRB1*07:01,HLA-DPA1*02:01/DPB1*01:01,HLA-DRB1*09:01,HLA-DPA1*02:01/DPB1*05:01,HLA-DRB3*01:01,HLA-DPA1*03:01/DPB1*04:02,HLA-DRB4*01:01,HLA-DQA1*01:01/DQB1*05:01,HLA-DQA1*01:02/DQB1*06:02,HLA-DQA1*03:01/DQB1*03:02,HLA-DQA1*04:01/DQB1*04:02,HLA-DQA1*05:01/DQB1*02:01,HLA-DQA1*05:01/DQB1*03:01'
allele_consensus3='HLA-DPA1*01/DPB1*04:01,HLA-DRB1*01:01,HLA-DRB1*04:04,HLA-DRB1*08:04,HLA-DRB1*11:21,HLA-DRB1*13:27,HLA-DPA1*01:03/DPB1*02:01,HLA-DRB1*01:02,HLA-DRB1*04:05,HLA-DRB1*08:06,HLA-DRB1*11:28,HLA-DRB1*13:28,HLA-DPA1*02:01/DPB1*01:01,HLA-DRB1*03:01,HLA-DRB1*04:08,HLA-DRB1*08:13,HLA-DRB1*13:01,HLA-DRB1*15:01,HLA-DPA1*02:01/DPB1*05:01,HLA-DRB1*03:05,HLA-DRB1*04:10,HLA-DRB1*08:17,HLA-DRB1*13:02,HLA-DRB1*15:02,HLA-DPA1*03:01/DPB1*04:02,HLA-DRB1*03:06,HLA-DRB1*04:21,HLA-DRB1*11:01,HLA-DRB1*13:04,HLA-DRB1*15:06,HLA-DQA1*01:01/DQB1*05:01,HLA-DRB1*03:07,HLA-DRB1*04:23,HLA-DRB1*11:02,HLA-DRB1*13:05,HLA-DRB3*01:01,HLA-DQA1*01:02/DQB1*06:02,HLA-DRB1*03:08,HLA-DRB1*04:26,HLA-DRB1*11:04,HLA-DRB1*13:07,HLA-DRB4*01:01,HLA-DQA1*03:01/DQB1*03:02,HLA-DRB1*03:09,HLA-DRB1*07:01,HLA-DRB1*11:06,HLA-DRB1*13:11,HLA-DRB5*01:01,HLA-DQA1*04:01/DQB1*04:02,HLA-DRB1*03:11,HLA-DRB1*07:03,HLA-DRB1*11:07,HLA-DRB1*13:21,HLA-DRB5*01:05,HLA-DQA1*05:01/DQB1*02:01,HLA-DRB1*04:01,HLA-DRB1*08:01,HLA-DRB1*11:14,HLA-DRB1*13:22,HLA-DQA1*05:01/DQB1*03:01,HLA-DRB1*04:02,HLA-DRB1*08:02,HLA-DRB1*11:20,HLA-DRB1*13:23'

def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (description='IEDB classII prediction')
	parser.add_argument ('-o', dest='out_dir', help='The output path')
	parser.add_argument ('-a', dest='hla_file', help='hla_file')
	parser.add_argument ('-p', dest='pep_fasta', help='Input peptide fasta')
	parser.add_argument ('-f', dest='sample_flag', help='hla flag')
	args = parser.parse_args ()
	return args

def get_NetMHCIIpan_EL(length,out_dir1,pep_fasta,HLA_file,sample_flag):
	fin = open (HLA_file)
	alleles = fin.readline()
	hla_allele = alleles.strip ("\n").split (",")
	hlas=allele_netmhcII.split(",")
	cmd = ''
	for hla in hla_allele:
		if ("/") in hla:
			hla1_tmp1 = str (hla).split ("/")[0] + str (hla).split ("/")[1]
			hla1_tmp = str (hla1_tmp1).split ("*")[0] + str (hla1_tmp1).split ("*")[1] + str (hla1_tmp1).split ("*")[2]
			if (hla1_tmp.count (':') == 2):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1] + str (hla1_tmp).split (":")[2]
			elif (hla1_tmp.count (':') == 1):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		else:
			print(hla)
			if hla=='':
				continue
			else:
				hla1_tmp = str (hla).split ("*")[0] + str (hla).split ("*")[1]
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		cmd += IEDB + ' netmhciipan_el ' + str (hla) + ' ' + pep_fasta + ' ' + str(length) + ' '
		cmd += ' > ' + out_dir1 + '/' + sample_flag + '_IEDB_NetMHCIIpan_el_' + str (hla1) + '_length_' + str(length) + '.log' + '\n'
	return (cmd)

def get_NetMHCIIpan_BA(length,out_dir1,pep_fasta,HLA_file,sample_flag):
	fin = open (HLA_file)
	alleles = fin.readline()
	hla_allele = alleles.strip ("\n").split (",")
	hlas=allele_netmhcII.split(",")
	cmd = ''
	for hla in hla_allele:
		if ("/") in hla:
			hla1_tmp1 = str (hla).split ("/")[0] + str (hla).split ("/")[1]
			hla1_tmp = str (hla1_tmp1).split ("*")[0] + str (hla1_tmp1).split ("*")[1] + str (hla1_tmp1).split ("*")[2]
			if (hla1_tmp.count (':') == 2):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1] + str (hla1_tmp).split (":")[2]
			elif (hla1_tmp.count (':') == 1):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		else:
			print (hla)
			if hla == '':
				continue
			else:
				hla1_tmp = str (hla).split ("*")[0] + str (hla).split ("*")[1]
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		cmd += IEDB + ' netmhciipan_ba ' + str (hla) + ' ' + pep_fasta + ' ' + str(length) + ' '
		cmd += ' > ' + out_dir1 + '/' + sample_flag + '_IEDB_NetMHCIIpan_ba_' + str (hla1) + '_length_' + str(length) + '.log' + '\n'
	return (cmd)

def get_sturniolo(length,out_dir1,pep_fasta,HLA_file ,sample_flag):
	fin=open(HLA_file)
	alleles=fin.readline()
	hla_allele=alleles.strip("\n").split(",")
	hlas=allele_sturniolo.split(",")
	cmd = ''
	for hla in hla_allele:
		if ("/") in hla:
			continue
		else:
			print (hla)
			if hla == '':
				continue
			else:
				if hla in hlas:
					hla1_tmp = str (hla).split ("*")[0] + str (hla).split ("*")[1]
					hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
				cmd += IEDB + ' sturniolo ' + str (hla) + ' ' + pep_fasta + ' ' + str (length) + ' '
				cmd += ' > ' + out_dir1 + '/' + sample_flag + '_IEDB_sturniolo_' + str (hla1) + '_length_' + str (length) + '.log' + '\n'
	return (cmd)

def get_NetMHCIIpan(length,out_dir1,pep_fasta,HLA_file,sample_flag):
	fin = open (HLA_file)
	alleles = fin.readline()
	hla_allele = alleles.strip ("\n").split (",")
	hlas=allele_netmhcII.split(",")
	cmd = ''
	for hla in hla_allele:
		if ("/") in hla:
			hla1_tmp1 = str (hla).split ("/")[0] + str (hla).split ("/")[1]
			hla1_tmp = str (hla1_tmp1).split ("*")[0] + str (hla1_tmp1).split ("*")[1] + str (hla1_tmp1).split ("*")[2]
			if (hla1_tmp.count (':') == 2):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1] + str (hla1_tmp).split (":")[2]
			elif (hla1_tmp.count (':') == 1):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		else:
			print (hla)
			if hla == '':
				continue
			else:
				hla1_tmp = str (hla).split ("*")[0] + str (hla).split ("*")[1]
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		cmd += IEDB + ' NetMHCIIpan ' + str (hla) + ' ' + pep_fasta + ' ' + str(length) + ' '
		cmd += ' > ' + out_dir1 + '/' + sample_flag + '_IEDB_NetMHCIIpan_' + str (hla1) + '_length_' + str(length) + '.log' + '\n'
	return (cmd)


def get_nn(length,out_dir1,pep_fasta,HLA_file,sample_flag):
	fin = open (HLA_file)
	alleles = fin.readline()
	hla_allele = alleles.strip ("\n").split (",")
	hlas = allele_nn_align.split (",")
	cmd = ''
	for hla in hla_allele:
		if ("/") in hla:
			hla1_tmp1 = str (hla).split ("/")[0] + str (hla).split ("/")[1]
			hla1_tmp = str (hla1_tmp1).split ("*")[0] + str (hla1_tmp1).split ("*")[1] + str (hla1_tmp1).split ("*")[2]
			if (hla1_tmp.count (':') == 2):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1] + str (hla1_tmp).split (":")[2]
			elif (hla1_tmp.count (':') == 1):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		else:
			print (hla)
			if hla == '':
				continue
			else:
				hla1_tmp = str (hla).split ("*")[0] + str (hla).split ("*")[1]
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		cmd += IEDB + ' nn_align ' + str (hla) + ' ' + pep_fasta + ' ' + str(length) + ' '
		cmd += ' > ' + out_dir1 + '/' + sample_flag + '_IEDB_nn_' + str (hla1) + '_length_' + str(length) + '.log' + '\n'
	return (cmd)

def get_smm(length,out_dir1,pep_fasta,HLA_file,sample_flag):
	fin = open (HLA_file)
	alleles = fin.readline()
	hla_allele = alleles.strip ("\n").split (",")
	hlas = allele_smm_align.split (",")
	cmd=''
	for hla in hla_allele:
		if ("/") in hla:
			hla1_tmp1 = str (hla).split ("/")[0] + str (hla).split ("/")[1]
			hla1_tmp = str (hla1_tmp1).split ("*")[0] + str (hla1_tmp1).split ("*")[1] + str (hla1_tmp1).split ("*")[2]
			if (hla1_tmp.count (':') == 2):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1] + str (hla1_tmp).split (":")[2]
			elif (hla1_tmp.count (':') == 1):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		else:
			print (hla)
			if hla == '':
				continue
			else:
				hla1_tmp = str (hla).split ("*")[0] + str (hla).split ("*")[1]
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		cmd += IEDB + ' smm_align ' + str (hla) + ' ' + pep_fasta + ' ' + str(length) + ' '
		cmd += ' > ' + out_dir1 + '/' + sample_flag + '_IEDB_smm_' + str (hla1) + '_length_' + str(length) + '.log' + '\n'
	return (cmd)


def get_comblib(length,out_dir1,pep_fasta,HLA_file,sample_flag):
	fin = open (HLA_file)
	alleles = fin.readline()
	hla_allele = alleles.strip ("\n").split (",")
	hlas= allele_comblib.split(",")
	cmd = ''
	for hla in hla_allele:
		if ("/") in hla:
			hla1_tmp1 = str (hla).split ("/")[0] + str (hla).split ("/")[1]
			hla1_tmp = str (hla1_tmp1).split ("*")[0] + str (hla1_tmp1).split ("*")[1] + str (hla1_tmp1).split ("*")[2]
			if (hla1_tmp.count (':') == 2):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1] + str (hla1_tmp).split (":")[2]
			elif (hla1_tmp.count (':') == 1):
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		else:
			print (hla)
			if hla == '':
				continue
			else:
				hla1_tmp = str (hla).split ("*")[0] + str (hla).split ("*")[1]
				hla1 = str (hla1_tmp).split (":")[0] + str (hla1_tmp).split (":")[1]
		cmd += IEDB + ' comblib ' + str (hla) + ' ' + pep_fasta + ' ' + str(length) + ' '
		cmd += ' > ' + out_dir1 + '/' + sample_flag + '_IEDB_comblib_' + str (hla1) + '_length_' + str(length) + '.log' + '\n'
	return (cmd)



def main():
	args = parse_arguments ()
	out_dir = args.out_dir
	HLA_file = args.hla_file
	sample_flag=args.sample_flag
	pep_fasta = args.pep_fasta
	sh_dir = out_dir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
	log_dir = out_dir + "/log/"  # The script log/alignment_start_end.sh will be executed

	if os.path.isdir (sh_dir) and os.path.isdir (log_dir):
		pass
	else:
		subprocess.call ('mkdir -p ' + out_dir + '/{sh,log}/', shell=True)

	fout = sh_dir + sample_flag + '_iedbII.sh'
	sh = open (fout, 'w')
	log = log_dir + sample_flag + '_iedbII.log'

	#subprocess.call ('mkdir -p ' + out_dir + '/' + sample_flag +'/', shell=True)
	out_dir1=out_dir

	cmd = ''
	cmd += 'echo starting running IEDB methods at `date` \n\n'
	for length in range(15,26):
		cmd += get_NetMHCIIpan_EL(length,out_dir1,pep_fasta,HLA_file,sample_flag)
		cmd += get_NetMHCIIpan_BA(length,out_dir1,pep_fasta,HLA_file,sample_flag)
		cmd += get_nn(length,out_dir1,pep_fasta,HLA_file ,sample_flag)
		cmd += get_smm(length,out_dir1,pep_fasta,HLA_file ,sample_flag)
		cmd += get_sturniolo(length,out_dir1,pep_fasta,HLA_file,sample_flag)
		cmd += get_comblib(length,out_dir1,pep_fasta,HLA_file,sample_flag)
		sh.write (cmd)
	sh.close ()
	cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
	subprocess.call (cmd, shell=True)

if __name__ == '__main__':
	main()

