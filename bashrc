# .bashrc

I_DIR=$HOME/works/bin


# Change according to your openfoam installation directory
#source ~/OpenFOAM/OpenFOAM-v1606+/etc/bashrc 
source /opt/OF/OpenFOAM-v1906/etc/bashrc

export I_DIR
export myBinDir=$I_DIR/bin

export PATH=$PATH:$I_DIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$I_DIR/bin



#~ export PATH=$PATH:$HOME/works/apps/macros
#~ export PATH=$PATH:$HOME/works/apps/macros/SP
#~ export PATH=$PATH:$HOME/works/apps/macros/preprocessing


