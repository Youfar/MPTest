#!/usr/bin/env python
# -*- coding: utf-8 -*-

# this program is using EM algorithm to estimate parameters of a mixture of bisulfite reads from two component
import random

READ = [[0 for col in range(2)] for row in range(9)]

## READ means a set of reads of an interval
## first column means the index of the short read (0-8)
## second column means the CpG site methylation status, 0 means unmethylated, 1 means methylated
## the value means the number of CpG site
## for example, READ[1][1] = 9 means READ 1 have 9 methylated CpG site

READ[0][0] = 9
READ[0][1] = 0
READ[1][0] = 0
READ[1][1] = 9
READ[2][0] = 1
READ[2][1] = 6
READ[3][0] = 0
READ[3][1] = 11
READ[4][0] = 1
READ[4][1] = 5
READ[5][0] = 8
READ[5][1] = 1
READ[6][0] = 1
READ[6][1] = 8
READ[7][0] = 8
READ[7][1] = 1
READ[8][0] = 0
READ[8][1] = 7
# there are 3 unmethylated short reads and 6 methylated short reads in this interval
# so the expected estimate minor ratio result of EM algorithm should be about 30%

def EM(READ):
    """
        EM algorithm. input is a set of reads of an interval, output is m1 and m2 and alpha
        m1 is methylation level of the minor part
        m2 is methylation level of the major part
        alpha is the minor part composition ratio
    """

    m1 = random.random()
    m2 = random.random() ## initial values of m1 and m2
    alpha1 = random.random()

    steps = 0 ## number of EM steps
    while True:
        nu1 = 0
        de1 = 0
        nu2 = 0
        de2 = 0

        nu_alpha = 0
        for key in range(9):
            ## one short read should have at least one CpG site
            if READ[key][1] + READ[key][0] >= 1:
                ## according to the formula
                p1 = (m1**READ[key][1])*((1-m1)**(READ[key][0]))
                p2 = (m2**READ[key][1])*((1-m2)**(READ[key][0]))
                Q1 = alpha1*p1 / (alpha1*p1 + (1-alpha1)*p2)
                Q2 = (1-alpha1)*p2 / (alpha1*p1 + (1-alpha1)*p2)

                ## nu1,nu2,de1,de2 is using for calculating m1,m2
                nu1 += Q1 * READ[key][1]
                de1 += Q1 * (READ[key][1] + READ[key][0])
                nu2 += Q2 * READ[key][1]
                de2 += Q2 * (READ[key][1] + READ[key][0])

                nu_alpha += Q1

        if (de1*de2 ==0): ## invalid, return a random result
            return (m1,m2,alpha1)
            break

        m1_new = nu1 / de1
        m2_new = nu2 / de2

        alpha1_new = nu_alpha/len(READ)

        if (abs(m1_new - m1) < 0.01 and abs(m2_new - m2) < 0.01 and abs(alpha1_new - alpha1) < 0.01) or steps > 200: ## num of steps <= 200
                if (abs(m1_new - m1) < 0.01 and abs(m2_new - m2) < 0.01 and abs(alpha1_new - alpha1) < 0.01):

                    # alpha1 is the smaller one
                    if (alpha1 > 0.5): # switch
                        alpha1 = 1 - alpha1
                        temp = m1
                        m1 = m2
                        m2 = temp
                    return (m1,m2,alpha1)
                    break
                elif steps > 200:
                    # alpha1 is the smaller one
                    if (alpha1 > 0.5): # switch
                        alpha1 = 1 - alpha1
                        temp = m1
                        m1 = m2
                        m2 = temp
                    return (m1,m2,alpha1)
                    break
        m1 = m1_new
        m2 = m2_new
        alpha1 = alpha1_new
        steps = steps + 1


def EM(READ):
    """ EM algorithm. input is a set of reads of an interval, output is m1 and m2 and alpha """
    #这个步骤是EM的第一次使用，主要用来求准确的alpha
    #这个设定说明结果和m1, m2, alpha初始值关系不大
    #这个算法根本没考虑，甲基化或者非甲基化的CpG site到底在short read的哪个位置
    #原来的实现是字典，
    m1 = random.random()
    m2 = random.random() ## initial values of m1 and m2
    alpha1 = random.random()

    steps = 0 ## number of EM steps
    while True:
        nu1 = 0
        de1 = 0
        nu2 = 0
        de2 = 0

        nu_alpha = 0
        for key in range(9):
            #一条短序列里面至少有一个CpG site
            if READ[key][1] + READ[key][0] >= 1:
                p1 = (m1**READ[key][1])*((1-m1)**(READ[key][0]))
                p2 = (m2**READ[key][1])*((1-m2)**(READ[key][0]))
                Q1 = alpha1*p1 / (alpha1*p1 + (1-alpha1)*p2)
                Q2 = (1-alpha1)*p2 / (alpha1*p1 + (1-alpha1)*p2)
                nu1 += Q1 * READ[key][1]
                de1 += Q1 * (READ[key][1] + READ[key][0])
                nu2 += Q2 * READ[key][1]
                de2 += Q2 * (READ[key][1] + READ[key][0])

                nu_alpha += Q1



        #print de1
        #print de2
        #sys.exit(0)
        if (de1*de2 ==0): ## invalid, return a random result
            return (m1,m2,alpha1)
            break

        m1_new = nu1 / de1
        m2_new = nu2 / de2

        alpha1_new = nu_alpha/len(READ)

        if (abs(m1_new - m1) < 0.01 and abs(m2_new - m2) < 0.01 and abs(alpha1_new - alpha1) < 0.01) or steps > 200: ## num of steps <= 200
                if (abs(m1_new - m1) < 0.01 and abs(m2_new - m2) < 0.01 and abs(alpha1_new - alpha1) < 0.01):

                    # alpha1 is the smaller one
                    if (alpha1 > 0.5): # switch
                        alpha1 = 1 - alpha1
                        temp = m1
                        m1 = m2
                        m2 = temp
                    return (m1,m2,alpha1)
                    break
                elif steps > 200:
                    #print "m1,m2,alpha1 not converge in 200 steps by EM calculating!\n"

                    # alpha1 is the smaller one
                    if (alpha1 > 0.5): # switch
                        alpha1 = 1 - alpha1
                        temp = m1
                        m1 = m2
                        m2 = temp
                    return (m1,m2,alpha1)
                    break
        m1 = m1_new
        m2 = m2_new
        alpha1 = alpha1_new
        steps = steps + 1

#def EM_with_alpha(READ, alpha1, m1_start, m2_start):

def bootstrap(READ, repeated_times):
    """
    其实就是用bootstrap反复求qualifying bins里面M1, M2, alpha1的方差
    可以弄出简易的结果
    :param READ:
    :param repeated_times:
    :return:
    """
    ##表示每个bins中的short read的数目
    num_of_read = len(READ)

    ## M1,M2,alpha1是每次循环参数的数组，用来求方差的
    M1 = [0]*repeated_times
    M2 = [0]*repeated_times
    alpha1 = [0]*repeated_times


    for i in range(repeated_times):
        ##SAMPLE_key表示从READ中选择short read的序号
        SAMPLE_key = []
        temp = []
        sample_times = min(1000, num_of_read)

        temp = [0,1,2,3,4,5,6,7,8]

        for j in range(sample_times):
            SAMPLE_key.append(random.randint(0, 8))

        #SAMPLE = {}
        SAMPLE = [[0 for col in range(2)] for row in range(9)]
        for k in SAMPLE_key:
            SAMPLE[k][0] = READ[k][0]
            SAMPLE[k][1] = READ[k][1]
        (m1_s, m2_s, al_s) = EM(SAMPLE)

        alpha1[i] = al_s
        M1[i] = m1_s
        M2[i] = m2_s
        del SAMPLE_key
        del SAMPLE

    var_alpha = numpy.var(alpha1)
    var_M1 = numpy.var(M1)
    var_M2 = numpy.var(M2)

    return (var_alpha, var_M1, var_M2)


m1,m2,alpha1 = EM(READ);
print m1, m2, alpha1



