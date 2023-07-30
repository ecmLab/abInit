#!/usr/bin/env python
######
#Script for parsing the DOS and visuliaze


__Author__='Bin Ouyang'
import os
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, Orbital
import numpy as np
from copy import deepcopy
import matplotlib as mpl
mpl.use('Agg');
from matplotlib import pyplot as plt
import argparse


def get_orbital_pdos(Energies,PDOS,NAtom,maxaqn=3,sm=2,ELim=[-5.0,5.0]):
    '''
    Get the orbital projected density of states,
    will return the figure handle for further process
    '''
    OrbDOSColor=['g','y','r','b'];
    NOrb=0; OrbIndsGrp=[]; Spins=[]; OrbLst=list(Orbital);
    OrbDict={0:'s',1:'p',2:'d',3:'f'};
    SpinName={0:Spin.up,1:Spin.down};
    for aqn in range(1,maxaqn+1):
        OrbNames=[OrbLst[i] for i in range(NOrb,NOrb+2*aqn-1)];
        OrbIndsGrp.append(OrbNames);
        NOrb+=2*aqn-1;
    for s in range(sm): Spins.append(SpinName[s]);

    ###Start the visualization
    PDOSArryDict={}; MaxDOS=0; MinDOS=1e10;
    Fig=plt.figure(figsize=(4,6));
    for GrpInd, OrbNames in enumerate(OrbIndsGrp):
        for si,s in enumerate(Spins):
            DOSArry=np.zeros(PDOS[0][OrbNames[0]][s].shape);
            for AtomInd in range(NAtom):
                for OrbName in OrbNames:
                    if si%2==0: DOSArry+=PDOS[AtomInd][OrbName][s];
                    else: DOSArry-=PDOS[AtomInd][OrbName][s];
            if si%2==0:
                plt.plot(DOSArry,Energies,color=OrbDOSColor[GrpInd],\
                    linewidth=2.0,label=OrbDict[GrpInd]);
                Tag='{}_up'.format(OrbDict[GrpInd]);
            else:
                plt.plot(DOSArry,Energies,color=OrbDOSColor[GrpInd],linewidth=2.0);
                Tag='{}_dn'.format(OrbDict[GrpInd]);
            PDOSArryDict[Tag]=DOSArry;
            if np.max(DOSArry[5:-5]) > MaxDOS: MaxDOS=np.max(DOSArry[5:-5]);
            if np.min(DOSArry[5:-5]) < MinDOS: MinDOS=np.min(DOSArry[5:-5]);
    ####Some logistics for visualization
    plt.xlabel('DOS (a.u.)',fontsize=21); plt.ylabel('E-E$_f$ (eV)',fontsize=21);
    plt.xlim((MinDOS*1.1,MaxDOS*1.1));
    plt.ylim((ELim[0],ELim[1]));
    plt.xticks(fontsize=18); plt.yticks(fontsize=18);
    plt.legend(fontsize=18,loc='best');
    plt.tight_layout();
    return plt, PDOSArryDict;

def get_atom_pdos(EleIndDict,EleSyms,Energies,PDOS,NAtom,maxsqn=3,sm=2,ELim=[-5.0,5.0]):
    '''
    Get the atomic projected density of states,
    will return the figure handle for further process
    '''
    AtomDOSColor = ['g','y','r','b'];
    NOrb=0; OrbIndsGrp=[]; Spins=[]; OrbLst=list(Orbital);
    OrbDict={0:'s',1:'p',2:'d',3:'f'};
    SpinName={0:Spin.up,1:Spin.down};
    for s in range(sm): Spins.append(SpinName[s]);

    ###Start the visualization
    PDOSArryDict={}; MaxDOS=0; MinDOS=1e10;
    Fig=plt.figure(figsize=(4,6));
    for EleInd, EleSym in enumerate(EleSyms):
        IndLst=EleIndDict[EleSym];
        for si,s in enumerate(Spins):
            DOSArry=np.zeros(PDOS[0][OrbLst[0]][s].shape);
            for AtomInd in IndLst:
                for OrbInd,Orb in enumerate(OrbLst):
                    if si%2==0: DOSArry+=PDOS[AtomInd][Orb][s];
                    else: DOSArry-=PDOS[AtomInd][Orb][s];
            if si%2==0:
                plt.plot(DOSArry,Energies,color=AtomDOSColor[EleInd],\
                        linewidth=2.0,label=EleSym);
                Tag='{}_up'.format(EleSym);
            else:
                plt.plot(DOSArry,Energies,color=AtomDOSColor[EleInd],linewidth=2.0);
                Tag='{}_dn'.format(EleSym);
            PDOSArryDict[Tag]=DOSArry;
            if np.max(DOSArry[5:-5]) > MaxDOS: MaxDOS=np.max(DOSArry[5:-5]);
            if np.min(DOSArry[5:-5]) < MinDOS: MinDOS=np.min(DOSArry[5:-5]);
    ####Some logistics for visualization
    plt.xlabel('DOS (a.u.)',fontsize=21); plt.ylabel('E-E$_f$ (eV)',fontsize=21);
    plt.xlim((MinDOS*1.1,MaxDOS*1.1));
    plt.ylim((ELim[0],ELim[1]));
    plt.xticks(fontsize=18); plt.yticks(fontsize=18);
    plt.legend(fontsize=18,loc='best');
    plt.tight_layout();
    return plt, PDOSArryDict;

def extractVasprun(VPRunFile):
    '''
    Extract neccessary information from vasprun.xml
    '''
    VPR=Vasprun(VPRunFile);
    TDOS,IDOS,PDOS=VPR.tdos,VPR.idos,VPR.pdos;
    Energies=TDOS.energies-TDOS.efermi;
    ProtoStr=VPR.structures[-1];
    EleSyms=list(ProtoStr.symbol_set);
    EleIndDict={Ele:list(ProtoStr.indices_from_symbol(Ele)) for Ele in EleSyms};
    NAtom=len(ProtoStr);

    return Energies, PDOS, NAtom, EleIndDict, EleSyms;


if __name__ == "__main__":
    Parser=argparse.ArgumentParser();
    Parser.add_argument('--vprun',type=str,help='Vasprun file (default: vasprun.xml)',\
            default='vasprun.xml');
    Parser.add_argument('--aqn',type=int,default=3,\
            help='Highest angular quantum number to be parsed (default: 3 (up to d))');
    Parser.add_argument('--sm',type=int,default=2,\
            help='Highest spin (default: 2 (spin polarized))');
    Parser.add_argument('--orb',help="Get orbital projected DOS",action='store_true');
    Parser.add_argument('--atom',help="Get atomic projected DOS",action='store_true');
    Parser.add_argument('--elim',help="The ylim/elim ",nargs='+',type=float);
    Parser.add_argument('--einv',help="The interval between two energy",type=float,default=2.0);
    Parser.add_argument('--figName',help="The figure name (default: PDOS.pdf)",type=str,\
            default='PDOS.pdf');
    Args=Parser.parse_args();

    Energies,PDOS,NAtom,EleIndDict,EleSyms=extractVasprun(Args.vprun);

    if Args.orb:
        plt,aaa=get_orbital_pdos(Energies,PDOS,NAtom,Args.aqn,Args.sm,Args.elim);
        plt.tight_layout(); plt.savefig(Args.figName);
        print(aaa);
        

    if Args.atom:
        plt,_=get_atom_pdos(EleIndDict,EleSyms,Energies,PDOS,NAtom,Args.aqn,Args.sm,Args.elim);
        plt.tight_layout(); plt.savefig(Args.figName);


