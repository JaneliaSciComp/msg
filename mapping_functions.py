##mapping functions
##0.2.1 Now passes barcodes_file name supplied by user.

import re
import os
import gzip
try:
    from pylab import *
    from matplotlib.font_manager import FontProperties
except:
    pass ## Dan: so that this doesn't create noise during main analysis run

def read_barcodes(barcodes_file):
    barcodes_file = open(barcodes_file,'r')
#    barcodes_file.readline()##ignore first line
#    variables = barcodes_file.readline().strip('\n').split('\t')
    bc = []
    for line in barcodes_file:
        y = line.strip('\n').split()
        #y = line.strip('\n').split('\t')
        bc.append(y) #y[n][0] is the barcode, y[n][1] is the indiv identifier
    barcodes_file.close()
    return bc
#    return variables, bc

def count_lines(count_file):
    if not os.path.exists(count_file):
        count_file += '.gz'
    try:
        if count_file.lower().endswith('.gz'):
            filename = gzip.open(count_file, 'rb')
        else:
            filename = open(count_file,'r')
    except IOError, Err:
        print "count lines of '%s'. %s" % (count_file, str(Err))
        return 0
    x = 0
    for line in filename:
        x+=1
    filename.close()
    return x

def create_stats(raw_data, barcodes_file):
    """Write out miscellaneous parsing stats."""
    # "indiv#_barcode"
    stats_file_name = raw_data+'_stats.txt'
    if os.path.exists(stats_file_name):
        print "Stats file '%s' already exists. Skipping" % stats_file_name
        return
    stats_file = open(stats_file_name,'w')
    #print file with statistics of parsing
    total_reads = 0
    read_file = './'+raw_data+'_parsed/bad_barcodes'
    if os.path.exists(read_file): #Catch the case where the bad_barcodes file wasn't created or copied over to the parsed reads directory
        number_reads = count_lines(read_file)/4
        stats_file.write('bad_barcodes\t%s\n' %(number_reads))
        total_reads += number_reads
    else:
        stats_file.write('bad_barcodes\tNA\n')

    read_file = './'+raw_data+'_parsed/unreadable_barcodes'
    if os.path.exists(read_file): #Catch the case where the unreadable_barcodes file wasn't created or copied over to the parsed reads directory
        number_reads = count_lines(read_file)/4
        stats_file.write('unreadable_barcodes\t%s\n' %(number_reads))   
        total_reads += number_reads
    else:
        stats_file.write('unreadable_barcodes\tNA\n')

    read_file = './'+raw_data+'_parsed/linkers'
    if os.path.exists(read_file): #Catch the case where the linkers file wasn't created or copied over to the parsed reads directory
        number_reads = count_lines(read_file)/2
        stats_file.write('linkers\t%s\n' %(number_reads))   
        total_reads += number_reads
    else:
        stats_file.write('linkers\tNA\n')

    read_file = './'+raw_data+'_parsed/junk'
    if os.path.exists(read_file): #Catch the case where the junk file wasn't created or copied over to the parsed reads directory
        number_reads = count_lines(read_file)/4
        stats_file.write('junk\t%s\n' %(number_reads))   
        total_reads += number_reads
    else:
        stats_file.write('junk\tNA\n')

    total_reads = 0 #Why?
    for ind in read_barcodes(barcodes_file):
        file_name = 'indiv' + ind[1] + '_' + ind[0]
        fastq_file = raw_data + '_parsed/' + file_name
        try:
            number_reads = count_lines(fastq_file)/4
        except IOError:
            #could have been gzipped
            number_reads = count_lines(fastq_file+'.gz')/4
        total_reads += number_reads
        stats_file.write('indiv%s_%s\t%s\n' %(ind[1],ind[0],number_reads))
   
    stats_file.write('total_reads\t%s' %(total_reads))
    stats_file.close()

def call_genotypes(data_filename,geno_span,ind_list,cross_type,sex,input_file,sld_wdw,sld_wdw_step,cutoffs,sp1,sp2,chrom_lgth):
    indiv_name = "indiv"+ ind_list[1] + "_" + ind_list[0]
    number_markers = []
    output_file = open(data_filename + "_genotypes_" + str(geno_span/1000) + 'kb/'+indiv_name + '_genotype.txt','w')
    data_file = open(input_file,'r')

    #create lists for genotype calls
    X_geno = []
    TwoL_geno = []
    TwoR_geno = []
    ThreeL_geno = []
    ThreeR_geno = []
    Four_geno = []

    #create lists for each chromosome arm
    X = []
    TwoL = []
    TwoR = []
    ThreeL = []
    ThreeR = []
    Four = []

    #read input file, parse chromosome data into 6 different arrays containing chromosome location & genotype call


    for line in data_file:
        x = line.strip('\n')
        data = x.split('\t')
        if data == ['']:#exits if extra empty lines are accidently pasted to the end of file
            pass
        elif "X" in data[2] and "chr" not in data[2]:
            X.append((int(data[3]),data[0]))
        elif "2L" in data[2] and "chr" not in data[2]:
            TwoL.append((int(data[3]),data[0]))
        elif "2R" in data[2] and "chr" not in data[2]:
            TwoR.append((int(data[3]),data[0]))
        elif "3L" in data[2] and "chr" not in data[2]:
            ThreeL.append((int(data[3]),data[0]))
        elif "3R" in data[2] and "chr" not in data[2]:
            ThreeR.append((int(data[3]),data[0]))
        elif "4" in data[2] and "chr" not in data[2]:
            Four.append((int(data[3]),data[0]))
        else:pass


##    #count up # markers
    number_markers.extend([len(X),len(TwoL),len(TwoR),len(ThreeL),len(ThreeR),len(Four),len(X+TwoL+TwoR+ThreeL+ThreeR+Four)])
##
##    #pull out unique genotype calls
    X = list(set(X))
    TwoL = list(set(TwoL))
    TwoR = list(set(TwoR))
    ThreeL = list(set(ThreeL))
    ThreeR = list(set(ThreeR))
    Four = list(set(Four))
    #sort list by chromosome location
    X.sort()
    TwoL.sort()
    TwoR.sort()
    ThreeL.sort()
    ThreeR.sort()
    Four.sort()
    #count up # markers after reduction
    number_markers.extend([len(X),len(TwoL),len(TwoR),len(ThreeL),len(ThreeR),len(Four),len(X+TwoL+TwoR+ThreeL+ThreeR+Four)])
    #get marker location and score
    #summarize over intervals of geno_span

    ##Print out file of unique sorted genotype calls
    ##Use this for future utility program that prints out portion of genome

    outfile = open("./"+data_filename+"_sorted_unique_genotypes/" + indiv_name + "_sorted_unique.txt","w")
    for line in X:
        outfile.write('%s\t%s\t%s\n' %("X",str(line[0]),str(line[1])))
    for line in TwoL:
        outfile.write('%s\t%s\t%s\n' %("2L",str(line[0]),str(line[1])))
    for line in TwoR:
        outfile.write('%s\t%s\t%s\n' %("2R",str(line[0]),str(line[1])))
    for line in ThreeL:
        outfile.write('%s\t%s\t%s\n' %("3L",str(line[0]),str(line[1])))
    for line in ThreeR:
        outfile.write('%s\t%s\t%s\n' %("3R",str(line[0]),str(line[1])))
    for line in Four:
        outfile.write('%s\t%s\t%s\n' %("4",str(line[0]),str(line[1])))
    outfile.close()
                      
    #get each marker
#    print X[0]
    X_geno = get_genotype_v2(X,geno_span,"X",cutoffs,cross_type,sex,sp1,sp2,chrom_lgth['X'])
    TwoL_geno = get_genotype_v2(TwoL,geno_span,"2L",cutoffs,cross_type,sex,sp1,sp2,chrom_lgth['2L'])
    TwoR_geno = get_genotype_v2(TwoR,geno_span,"2R",cutoffs,cross_type,sex,sp1,sp2,chrom_lgth['2R'])
    ThreeL_geno = get_genotype_v2(ThreeL,geno_span,"3L",cutoffs,cross_type,sex,sp1,sp2,chrom_lgth['3L'])
    ThreeR_geno = get_genotype_v2(ThreeR,geno_span,"3R",cutoffs,cross_type,sex,sp1,sp2,chrom_lgth['3R'])
    Four_geno = get_genotype_v2(Four,geno_span,"4",cutoffs,cross_type,sex,sp1,sp2,chrom_lgth['4'])

    ##print 'X\r',X_geno
    ##print '2L\r',TwoL_geno
    ##print '2R\r',TwoR_geno
    ##print '3L\r',ThreeL_geno
    ##print '3R\r',ThreeR_geno
    ##print '4\r',Four_geno


    #write data to text file
    #columns headings
    output_file.write('id,')
    for data in X_geno:
        output_file.write('X:%s-%s,' %(data[0],data[1]))
    for data in TwoL_geno:
        output_file.write('2L:%s-%s,' %(data[0],data[1]))
    for data in TwoR_geno:
        output_file.write('2R:%s-%s,' %(data[0],data[1]))
    for data in ThreeL_geno:
        output_file.write('3L:%s-%s,' %(data[0],data[1]))
    for data in ThreeR_geno:
        output_file.write('3R:%s-%s,' %(data[0],data[1]))
    n = 0
    for data in Four_geno:
        n +=1
        if n < len(Four_geno):
            output_file.write('4:%s-%s,' %(data[0],data[1]))
        else:
            output_file.write('4:%s-%s' %(data[0],data[1]))
    output_file.write('\n')
    #chromosome number
    output_file.write(',')
    for data in X_geno:
        output_file.write('X,')
    for data in TwoL_geno:
        output_file.write('2,')
    for data in TwoR_geno:
        output_file.write('2,')
    for data in ThreeL_geno:
        output_file.write('3,')
    for data in ThreeR_geno:
        output_file.write('3,')
    n = 0
    for data in Four_geno:
        n +=1
        if n < len(Four_geno):
            output_file.write('4,')
        else:
            output_file.write('4')            
    output_file.write('\n')
    #genotype call
    #first print individual id
    output_file.write('%s,' %(ind_list[1]))
    for data in X_geno:
        output_file.write('%s,' %(data[2]))
    for data in TwoL_geno:
        output_file.write('%s,' %(data[2]))
    for data in TwoR_geno:
        output_file.write('%s,' %(data[2]))
    for data in ThreeL_geno:
        output_file.write('%s,' %(data[2]))
    for data in ThreeR_geno:
        output_file.write('%s,' %(data[2]))
    n = 0
    for data in Four_geno:
        n +=1
        if n < len(Four_geno):
            output_file.write('%s,' %(data[2]))
        else:
            output_file.write('%s' %(data[2]))            
    output_file.write('\n')
    #sim calls
    output_file.write(sp1 + ' matches,')

    for data in X_geno:
        output_file.write('%s,' %(data[3]))
    for data in TwoL_geno:
        output_file.write('%s,' %(data[3]))
    for data in TwoR_geno:
        output_file.write('%s,' %(data[3]))
    for data in ThreeL_geno:
        output_file.write('%s,' %(data[3]))
    for data in ThreeR_geno:
        output_file.write('%s,' %(data[3]))
    n = 0
    for data in Four_geno:
        n +=1
        if n < len(Four_geno):
            output_file.write('%s,' %(data[3]))
        else:
            output_file.write('%s' %(data[3]))            
    output_file.write('\n')
    #sec calls
    output_file.write(sp2 + ' matches,')
    for data in X_geno:
        output_file.write('%s,' %(data[4]))
    for data in TwoL_geno:
        output_file.write('%s,' %(data[4]))
    for data in TwoR_geno:
        output_file.write('%s,' %(data[4]))
    for data in ThreeL_geno:
        output_file.write('%s,' %(data[4]))
    for data in ThreeR_geno:
        output_file.write('%s,' %(data[4]))
    n = 0
    for data in Four_geno:
        n +=1
        if n < len(Four_geno):
            output_file.write('%s,' %(data[4]))
        else:
            output_file.write('%s' %(data[4]))            


    output_file.close()
    data_file.close()


    #print out detailed maps of genotype calls 
    print_maps(ind_list,sld_wdw,sld_wdw_step,X,TwoL,TwoR,ThreeL,ThreeR,Four,sp1,sp2,chrom_lgth)
    #move map to plots folder
    genotype_plot = 'indiv' + ind_list[1] + '_' + ind_list[0] + "_plot.png"
    os.rename(genotype_plot, './' + data_filename + '_plots_' + str(sld_wdw/1000) + 'kb/' + genotype_plot)

    return number_markers




def get_data_in_interval(data,start,stop):
    x = []
    for z in data:
        if z[0] < start:
            pass
        elif z[0] >= start and z[0] < stop: #if datum in interval then add to list
            x.append(z[1])
        elif z[0] >= stop: #if we have gone beyond the interval 
            break            
    return x

def accum_data_in_interval(data,sp1,sp2,start,stop):
    sp1_count = 0
    sp2_count = 0
    x = []
    for z in data:
        if z[0] < start:
            pass
        elif z[0] >= start and z[0] < stop: #if datum in interval then add to list
            if sp1 in z[1]:sp1_count +=1
            elif sp2 in z[1]:sp2_count +=1
        elif z[0] >= stop: #if we have gone beyond the interval 
            break            


    return sp1_count,sp2_count


def get_genotype_v2(sorted_array_of_all_markers,step,chromosome,cutoffs,cross_type,sex,sp1,sp2,chrom_lgth):
    genotype_array=[]             
    start = 0
    stop = step
    num_steps = (chrom_lgth + step)/step #last interval should include the chromosome length, add step size, then divide by step to get integer # steps

    if sorted_array_of_all_markers: #check to make sure there is data in this array
        #check for number of markers in each interval, use this list to scan through data
        next_step = 0
        while next_step < num_steps:
            next_step = next_step + 1
            
            sp1_count,sp2_count = accum_data_in_interval(sorted_array_of_all_markers,sp1,sp2,start,stop)
                
            if sp1_count == 0 and sp2_count == 0:#if no data in interval
                genotype = '?'
                genotype_array.append([start, stop, genotype,sp1_count,sp2_count])
            else:#if some data in interval
                ratio = float(sp1_count-sp2_count)/(sp1_count+sp2_count)
                if cross_type == 'BCsp1':
                    if sex == 'male' and chromosome == 'X':
                        if ratio >= cutoffs['XBCmalesp1']:
                            genotype = sp1
                        else:
                            genotype = sp2
                    elif sex == 'female' and chromosome == 'X': 
                        if ratio >= cutoffs['XBCfemalesp1']:
                            genotype = sp1
                        else:
                            genotype = 'het'
                    elif chromosome != 'X':
                        if ratio >= cutoffs['AutBCsp1']:
                            genotype = sp1
                        else:
                            genotype = 'het'
                    else:#if any other symbol used, like ? or NA, then assume male for safer cutoff
                        if ratio >= cutoffs['XBCmalesp1']:
                            genotype = sp1
                        else:
                            genotype = sp2

                elif cross_type == 'BCsp2':
                    
                    if sex == 'male' and chromosome == 'X':
                        if ratio >= cutoffs['XBCmalesp2']:
                            genotype = sp1
                        else:
                            genotype = sp2
                    elif sex == 'female' and chromosome == 'X': 
                        if ratio >= cutoffs['XBCfemalesp2']:
                            genotype = 'het'
                        else:
                            genotype = sp2
                    elif chromosome != 'X':
                        if ratio >= cutoffs['AutBCsp2']:
                            genotype = 'het'
                        else:
                            genotype = sp2
                    else:#if any other symbol used, like ? or NA, then assume male for safer cutoff
                        if ratio >= cutoffs['XBCmalesp2']:
                            genotype = sp1
                        else:
                            genotype = sp2                            
                genotype_array.append([start, stop, genotype,sp1_count,sp2_count])

            stop = stop + step
            start = start + step

        return genotype_array
    
    else:#if no data in the array, fill with ???
        start = 0
        stop = step
        next_step = 0
        while next_step < num_steps:
            next_step = next_step + 1
            genotype_array.append([start,stop,'?',0,0])
            stop = stop + step
            start = start + step
        return genotype_array



##
##Plots over sliding window, not running average
##
def print_maps(ind_list,sld_wdw,sld_wdw_step,X,TwoL,TwoR,ThreeL,ThreeR,Four,sp1,sp2,chrom_lgth):
    #These lists contain [map position, genotype call]
    #want to replace genotype calls with numbers
#    print X[0:100]
    X_nums = map(lambda y: [y[0],1] if y==(y[0],sp1) else [y[0],-1], X)
    TwoL_nums = map(lambda y: [y[0],1] if y==(y[0],sp1) else [y[0],-1], TwoL)
    TwoR_nums = map(lambda y: [y[0],1] if y==(y[0],sp1) else [y[0],-1], TwoR)
    ThreeL_nums = map(lambda y: [y[0],1] if y==(y[0],sp1) else [y[0],-1], ThreeL)
    ThreeR_nums = map(lambda y: [y[0],1] if y==(y[0],sp1) else [y[0],-1], ThreeR)
    Four_nums = map(lambda y: [y[0],1] if y==(y[0],sp1) else [y[0],-1], Four)

    #plot each chromosome with running average
    #plot results

    subplots_adjust(wspace = 0.4, hspace = 0.4)
    suptitle('individual ' + ind_list[1] + '  barcode ' + ind_list[0], fontsize = 12)
    
    if X_nums:
        subplot(321)
        x=[]
        y=[]
        for datum in X_nums:
            x.append(datum[0])
            y.append(int(datum[1]))
        if len(x) > 2:#if more than 2 data points on chromosome
            run_points = sliding_window(X_nums,sld_wdw,sld_wdw_step,1,chrom_lgth['X'])
            plot(range(1+sld_wdw/2,chrom_lgth['X']-sld_wdw/2,sld_wdw_step),run_points,linewidth=1)
        scatter(x,y,s=10,marker=[2,0,0],c='black')
        axis([-1000,chrom_lgth['X']+1000,-1.1,1.1])
        xlabel('X')
       
    if Four_nums:
        subplot(322)
        x=[]
        y=[]
        for datum in Four_nums:
            x.append(datum[0])
            y.append(int(datum[1]))
        if len(x) > 2:#if more than 2 data points on chromosome
            run_points = sliding_window(Four_nums,sld_wdw,sld_wdw_step,1,chrom_lgth['4'])
            plot(range(1+sld_wdw/2,chrom_lgth['4']-sld_wdw/2,sld_wdw_step),run_points,linewidth=1)
        scatter(x,y,s=10,marker=[2,0,0],c='black')
        axis([-100,chrom_lgth['4']+100,-1.1,1.1])
        axhline(linewidth = 1, color = 'r')
        xlabel('4')

    if TwoL_nums:
        subplot(323)
        x=[]
        y=[]
        for datum in TwoL_nums:
            x.append(datum[0])
            y.append(int(datum[1]))
        if len(x) > 2:#if more than 2 data points on chromosome
            run_points = sliding_window(TwoL_nums,sld_wdw,sld_wdw_step,1,chrom_lgth['2L'])
#            print run_points
            plot(range(1+sld_wdw/2,chrom_lgth['2L']-sld_wdw/2,sld_wdw_step),run_points,linewidth=1)
        scatter(x,y,s=10,marker=[2,0,0],c='black')
        axhline(linewidth = 1, color = 'r')
        axis([-1000,chrom_lgth['2L']+1000,-1.1,1.1])
        xlabel('2L')

    if TwoR_nums:
        subplot(324)
        x=[]
        y=[]
        for datum in TwoR_nums:
            x.append(datum[0])
            y.append(int(datum[1]))
        if len(x) > 2:#if more than 2 data points on chromosome
            run_points = sliding_window(TwoR_nums,sld_wdw,sld_wdw_step,1,chrom_lgth['2R'])
            plot(range(1+sld_wdw/2,chrom_lgth['2R']-sld_wdw/2,sld_wdw_step),run_points,linewidth=1)
        scatter(x,y,s=10,marker=[2,0,0],c='black')
        axhline(linewidth = 1, color = 'r')
        axis([-1000,chrom_lgth['2R']+1000,-1.1,1.1])
        xlabel('2R')

    if ThreeL_nums:
        subplot(325)
        x=[]
        y=[]
        for datum in ThreeL_nums:
            x.append(datum[0])
            y.append(int(datum[1]))
        if len(x) > 2:#if more than 2 data points on chromosome
            run_points = sliding_window(ThreeL_nums,sld_wdw,sld_wdw_step,1,chrom_lgth['3L'])
            plot(range(1+sld_wdw/2,chrom_lgth['3L']-sld_wdw/2,sld_wdw_step),run_points,linewidth=1)
        scatter(x,y,s=10,marker=[2,0,0],c='black')
        axhline(linewidth = 1, color = 'r')
        axis([-1000,chrom_lgth['3L']+1000,-1.1,1.1])
        xlabel('3L')

    if ThreeR_nums:
        subplot(326)
        x=[]
        y=[]
        for datum in ThreeR_nums:
            x.append(datum[0])
            y.append(int(datum[1]))
        if len(x) > 2:#if more than 2 data points on chromosome
            run_points = sliding_window(ThreeR_nums,sld_wdw,sld_wdw_step,1,chrom_lgth['3R'])
            plot(range(1+sld_wdw/2,chrom_lgth['3R']-sld_wdw/2,sld_wdw_step),run_points,linewidth=1)
        scatter(x,y,s=10,marker=[2,0,0],c='black')
        axhline(linewidth = 1, color = 'r')
        axis([-1000,chrom_lgth['3R']+1000,-1.1,1.1])
        xlabel('3R')



    subplots_adjust(bottom=0.1, right = 0.8, top = 0.9)

    genotype_plot = 'indiv' + ind_list[1] + '_' + ind_list[0] + "_plot.png"
    savefig(genotype_plot, dpi = 300, format = 'png')
    close()
##    show()

##    output_file.close()


    
def sliding_window(data,sld_wdw,sld_wdw_step,start,end):#data is location,genotype
    num_steps = len(range(start+sld_wdw/2,end-sld_wdw/2,sld_wdw_step))
#    print num_steps
    next_step = 1
    stop = start + sld_wdw
    average = 0
    run_of_points = []
    while next_step <= num_steps:
        x = get_data_in_interval(data,start,stop)
        next_step = next_step + 1
        start = start + sld_wdw_step
        stop = stop + sld_wdw_step
        if x:
            average = float(sum(x))/len(x)
            run_of_points.append(average)#append average
        else:
            run_of_points.append(average)#if no data, append last average calculated from data
#    print len(run_of_points)
    return run_of_points
