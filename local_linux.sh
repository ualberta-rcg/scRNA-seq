WORKDIR=$PWD
mkdir lib
mkdir run
mkdir -p $WORKDIR/home/$USER
apptainer exec \
  --workdir $WORKDIR \
  --bind $WORKDIR/home/$USER:/home/$USER \
  --bind $WORKDIR/lib:/var/lib/rstudio-server \
  --bind $WORKDIR/run:/var/run/rstudio-server \
  --bind /etc/passwd \
  --bind /etc/group \
  --bind /tmp \
  --bind /var/tmp \
  scRNA-seq.sif \
  rserver --www-port=8787 --server-daemonize=0 --server-user=$USER
 