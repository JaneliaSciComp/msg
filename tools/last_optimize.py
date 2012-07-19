
"""
This was an attempt to optimize the parameters used for the Last alignment tool.
It didn't end up working. I believe I hit some overfitting issues.  To really run this right
you'd need it to run on several simulated files in the cost function.
Otherwise the code works though and could potentially be used to optimize paramaters for other
tools.

Greg Pinero 2012 (gregpinero@gmail.com)
"""

import mapping_benchmark

#Download the zip file from http://blog.kiwitobes.com/?p=44 (or search for programming collective 
#intelligence code) and copy optimization.py to this directory.  Change geneticoptimize function to
#return scores variable.
import optimization

def get_free_memory():
    """
    Try to figure out how much memory is free on a Unix system.
    Returns free memory in mB.
    """ 
    data = open("/proc/meminfo", 'rt').readlines()
    free = 0
    for line in data:
        if line.startswith("MemFree") or line.startswith("Buffers") or line.startswith("Cached"):
            items = line.split()
            free += int(items[1])
    #print "free memory is",free/1000,"mB"
    return free/1000

def gen_custom_args(sol):
    retval =''
    for i,(_,arg) in enumerate(DOMAIN):
        retval += ' -%s %s' % (arg, sol[i])
    return retval    

def cost(sol):
    """Lower cost is better"""
    #print "evaluating", gen_custom_args(sol)
    map = mapping_benchmark.LastForOptimize(mapping_benchmark.REF_FILE, 
        mapping_benchmark.SIM_READS_FILE)
    custom_args = gen_custom_args(sol)
    map.do_map(custom_args)
    map.filter_sam_file()
    run_time = map.get_run_time()
    #print "    run time: ", run_time
    #num_mapped = map.count_num_mapped()
    num_correct = map.calc_accuracy(ANSWERS)
    print "%s yields %s (free memory %s)" % (gen_custom_args(sol), num_correct, get_free_memory())
    return len(ANSWERS) - num_correct

def main():
    global ANSWERS
    global DOMAIN
    ANSWERS = mapping_benchmark.get_answers_by_qname(mapping_benchmark.SIM_READS_SAM_FILE)
    #define min and max values, and arg for the numeric arguments we want to optimize:
    DOMAIN = [
        ((1,300),'r'), #-r Score 
        ((1,300),'q'), #-q mismatch Cost
        ((1,300),'a'), #-a gap existence cost
        ((1,300),'b'), #-b gap extension cost
        #((1,300),'A'), #-A insertion existence cost, last says this is invalid??
        #((1,300),'B'), #-B insertion extension cost, last says this is invalid??
        #((),'c'), #-c gen affine gap costs? skip for now
        ((1,300),'x'), #-x drop
        ((1,300),'y'), #-y drop gaplees
        ((1,300),'z'), #-z drop final gapped
        ((1,300),'d'), #-d min score gapless
        ((1,300),'e'), #-e min alignment score
        ((1,300),'m'), #-m multipiclity
        ((30,100),'l'), #-l min length (less than 30 is way slow!)
        #((),'r'), #-n count skip for now?
        #((),'r'), #-t temperature?
        #((),'r'), #-g gamma?
    ]
    domain = [item[0] for item in DOMAIN]
    final_scores = optimization.geneticoptimize(domain=domain, costf=cost,
                                       popsize=300, step=3, mutprob=0.2, elite=0.2,
                                       maxiter=10)
    for (score,sol) in final_scores:
        print score, gen_custom_args(sol)

if __name__=='__main__':
    main()
