import os


def runqe(app,inf="qe.in",outf="qe.out"):
    """
    Run a quantum espresso calculation
    """
    from .calculations import cores,qepath
    if cores>1: os.system("mpirun -np "+str(cores)+" "+qepath+
            "/"+app+" -inp "+inf
            +" > "+outf)
    else: os.system(qepath+"/"+app+" < "+inf+" > "+outf)





def runelk(inf="elk.in",outf="elk.info"):
    """
    Run a quantum espresso calculation
    """
    from .calculations import cores,elkpath
    if cores>1: os.system("mpirun -np "+str(cores)+" "+elkpath+
            "/elk -inp "+inf
            +" > "+outf)
    else: os.system(elkpath+"/elk < "+inf+" > "+outf)
