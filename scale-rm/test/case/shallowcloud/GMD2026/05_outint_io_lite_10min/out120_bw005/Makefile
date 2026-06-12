################################################################################
# Makefile for GMD2026 shallow-cloud benchmark cases
################################################################################
PWD          = $(shell pwd)
TOPDIR       := $(abspath ../../../../../../..)
SCALE_RM_DIR := $(abspath ../../../../../..)
TESTDIR      := $(SCALE_RM_DIR)/test
index      = 200x11x11
precision  = real8
dynamics   = heve
tracer     = sn14
sfcflx     = const
turbulence = smg
radiation  = dycoms2
INITCONF  = init.conf
RUNCONF   = run.conf
TPROC     = 64
CODE_DIR = $(PWD)/code
ORG_SRCS = mod_user.f90
include $(TESTDIR)/Makefile.common
export SCALE_DYCOMS2_RF02_SDM=T
