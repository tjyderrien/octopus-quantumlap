#!gnuplot

set terminal postscript eps enhanced color font 'Helvetica, 24'

Efield="9E9" #3E9
CEP="0.0"
OrientationAngle="0"
kgrid1="24"
kgrid2="32"
kgrid3="40"
kgrid4="48"

set output 'HHGconvergence'.Efield.'Vm.eps'

filehhg="HHG-Spectrum.csv"
path1="/home/thibault.derrien/Documents/octopus.data/03-Bulk3D/Si-1x1x1-k-24x24x1-Angle0-Mesh0.434637008855-Mesh0.434637008855-T0e0/lambda2656e-9-tau30e-15-field".Efield."-CEP0.0-theta0-SinSquared-SinPower2-2/lambda400e-9-tau30e-15-field0E0-CEP0.0-delay0e-15-SinSquared/td.general/"
path2="/home/thibault.derrien/Documents/octopus.data/03-Bulk3D/Si-1x1x1-k-32x32x1-Angle0-Mesh0.434637008855-Mesh0.434637008855-T0e0/lambda2656e-9-tau30e-15-field".Efield."-CEP0.0-theta0-SinSquared-SinPower2-2/lambda400e-9-tau30e-15-field0E0-CEP0.0-delay0e-15-SinSquared/td.general/"
path3="/home/thibault.derrien/Documents/octopus.data/03-Bulk3D/Si-1x1x1-k-40x40x1-Angle0-Mesh0.4346370088553005-Mesh0.4346370088553005-T0e0/lambda2656e-9-tau13.945e-15-field".Efield."-CEP0.0-theta0-SinSquared-SinPower2-2/lambda400e-9-tau30e-15-field0E0-CEP0.0-delay0e-15-SinSquared/td.general/"
path4="/home/thibault.derrien/Documents/octopus.data/03-Bulk3D/Si-1x1x1-k-48x48x1-Angle0-Mesh0.4346370088553005-Mesh0.4346370088553005-T0e0/lambda2656e-9-tau13.945e-15-field".Efield."-CEP0.0-theta0-SinSquared-SinPower2-2/lambda400e-9-tau30e-15-field0E0-CEP0.0-delay0e-15-SinSquared/td.general/"

set xrange [:10]
set log y

set xlabel 'Energy (eV)'

p path1.filehhg u 1:2 w l t 'k=24', \
  path2.filehhg u 1:2 w l t 'k=32', \
  path3.filehhg u 1:2 w l t 'k=40', \
  path4.filehhg u 1:2 w l t 'k=48'
