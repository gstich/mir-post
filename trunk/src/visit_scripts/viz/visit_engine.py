# [VisIt Python Script]
# file: visit_engine.py
# author: Cyrus Harrison (cyrush@llnl.gov)
# created: 10/09/08
#
#
# Provides an Engine class that can be used to easily launch a visit
# engine on LLNL OCF machines.
#
# modifications:
#  Cyrus Harrison, Wed Jun 24 09:28:59 PDT 2009
#  Added hera.
#
#  Cyrus Harrison, Thu Jul  2 16:21:37 PDT 2009
#  Added prism.
#
#  Cyrus Harrison, Fri Jul 17 09:17:17 PDT 2009
#  Added yana.
#
 
import sys
import socket
import re
 
def parse_engine_args(args):
    """
    Helps parse optional engine args passed to the command line.
    """
    nargs  = len(args)
    nprocs = 8
    bank   = "<debug>"
    rtime  = "60:00"
    if nargs > 0:
        nprocs = int(args[0])
    if nargs > 1:
        bank = args[1]
    if nargs > 2:
        rtime = args[2]
    return nprocs, bank, rtime 
 
 
class Engine:
    def __init__(self, vdir = "/usr/gapps/visit/"):
        self.vdir = vdir
        self.host_name = self.__get_host_name()
    def open(self, nprocs = None, part = None, bank = None, rtime = None):
        lmethod_name = "_Engine__open_" + self.host_name
        lmethod = getattr(self,lmethod_name)
        lmethod(nprocs,part,bank,rtime)
    def close(self):
        """
        Closes VisIt's Compute Engine
        """
        CloseComputeEngine()
    def __get_host_name(self):
        host_name = socket.gethostname()
        return re.compile("[a-zA-z]*").match(host_name).group()
    def __proc_defaults(self,args,defs):
        nprocs = defs[0]
        part   = defs[1]
        bank   = defs[2]
        rtime  = defs[3]
        if args[0] is not None:
            nprocs = int(args[0])
        if args[1] is not None:
            part  = args[1]
        if args[2] is not None:
            bank  = args[2]
        if args[3] is not None:
            rtime = args[3]
        return nprocs, part, bank, rtime
    def __open_atlas(self,nprocs,part,bank,rtime):
        args = [nprocs,part,bank,rtime]
        defaults =[8,"pdebug","bdivp","60:00"]
        nprocs,part,bank,rtime = self.__proc_defaults(args, defaults)
        nnodes = nprocs / 8
        nnodes_str = "%s" % nnodes
        nprocs_str = "%s" % nprocs
        if part == "pbatch":
            print "[atlas: opening engine on pbatch]"
            OpenComputeEngine("localhost",("-np",nprocs_str,
                                           "-nn",nnodes_str,
                                           "-l", "msub/srun",
                                           "-b", bank,
                                           "-t", rtime,
                                           "-dir",self.vdir))
        else:
            print "[atlas: opening engine on pdebug]"
            OpenComputeEngine("localhost",("-np",nprocs_str,
                                           "-nn",nnodes_str,
                                           "-l", "srun",
                                           "-p", "pdebug",
                                           "-dir",self.vdir))
    def __open_hera(self,nprocs,part,bank,rtime):
        args = [nprocs,part,bank,rtime]
        defaults =[16,"pdebug","bdivp","60:00"]
        nprocs,part,bank,rtime = self.__proc_defaults(args, defaults)
        nnodes = nprocs / 16
        nnodes_str = "%s" % nnodes
        nprocs_str = "%s" % nprocs
        if part == "pbatch":
            print "[hera: opening engine on pbatch]"
            OpenComputeEngine("localhost",("-np",nprocs_str,
                                           "-nn",nnodes_str,
                                           "-l", "msub/srun",
                                           "-b", bank,
                                           "-t", rtime,
                                           "-dir",self.vdir))
        else:
            print "[hera: opening engine on pdebug]"
            OpenComputeEngine("localhost",("-np",nprocs_str,
                                           "-nn",nnodes_str,
                                           "-l", "srun",
                                           "-p", "pdebug",
                                           "-dir",self.vdir))
    def __open_yana(self,nprocs,part,bank,rtime):
        args = [nprocs,part,bank,rtime]
        defaults =[8,"pbatch","bdivp","30:00"]
        nprocs,part,bank,rtime = self.__proc_defaults(args, defaults)
        nnodes = nprocs / 8
        nnodes_str = "%s" % nnodes
        nprocs_str = "%s" % nprocs
        if part == "pbatch":
            print "[yana: opening engine on pbatch]"
            OpenComputeEngine("localhost",("-np",nprocs_str,
                                           "-nn",nnodes_str,
                                           "-l", "msub/srun",
                                           "-b", bank,
                                           "-t", rtime,
                                           "-dir",self.vdir))
        else:
            print "[yana: opening engine on current node]"
            OpenComputeEngine("localhost",("-np",nprocs_str,
                                           "-nn",nnodes_str,
                                           "-l", "srun",
                                           "-p", "pdebug",
                                           "-dir",self.vdir))
    def __open_zeus(self,nprocs,part,bank,rtime):
        args = [nprocs,part,bank,rtime]
        defaults =[8,"pdebug","bdivp","60:00"]
        nprocs,part,bank,rtime = self.__proc_defaults(args, defaults)
        nnodes = nprocs / 8
        nnodes_str = "%s" % nnodes
        nprocs_str = "%s" % nprocs
        if part == "pbatch":
            print "[zeus: opening engine on pbatch]"
            OpenComputeEngine("localhost",("-np",nprocs_str,
                                           "-nn",nnodes_str,
                                           "-l", "msub/srun",
                                           "-b", bank,
                                           "-t", rtime,
                                           "-dir",self.vdir))
        else:
            print "[zeus: opening engine on pdebug]"
            OpenComputeEngine("localhost",("-np",nprocs_str,
                                           "-nn",nnodes_str,
                                           "-l", "srun",
                                           "-p", "pdebug",
                                           "-dir",self.vdir))
    def __open_prism(self,nprocs,part,bank,rtime):
        args = [nprocs,part,bank,rtime]
        defaults =[8,"pbatch","bdivp","30:00"]
        nprocs,part,bank,rtime = self.__proc_defaults(args, defaults)
        nnodes = nprocs / 2
        nnodes_str = "%s" % nnodes
        nprocs_str = "%s" % nprocs
        if part == "pbatch":
            print "[prism: opening engine on pbatch]"
            OpenComputeEngine("localhost",("-np",nprocs_str,
                                           "-nn",nnodes_str,
                                           "-l", "msub/srun",
                                           "-b", bank,
                                           "-t", rtime,
                                           "-dir",self.vdir))
