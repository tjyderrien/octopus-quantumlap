# CONFIGURING VIDEO ON SALOMON. https://docs.it4i.cz/software/viz/NICEDCVsoftware/

1. Configure SSH tunnels
ssh -i ~/salomon_key -TN -f user@vizserv1.salomon.it4i.cz -L 5901:localhost:5901 -L 7300:localhost:7300 -L 7301:localhost:7301 -L 7302:localhost:7302 -L 7303:localhost:7303 -L 7304:localhost:7304 -L 7305:localhost:7305

2. Run the job
qsub -I -q qviz -A OPEN-16-27 -l select=1:ncpus=4:host=vizserv1,walltime=10:00:00 -N Vizserver1

3. Start X server
vncserver :1 -geometry 2540x1380 -depth 24
