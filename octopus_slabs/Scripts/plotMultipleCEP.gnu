#!gnuplot

CEP0="./lambda2237e-9-tau33.55e-15-dt0.25-field1.4E9-CEP0.0-theta0-SinFourth-SinPower2-2/lambda745.61e-9-tau33.55e-15-field7E7-CEP0.0-delay0e-15-SinFourth/td.general/laser"
CEPm05="./lambda2237e-9-tau33.55e-15-dt0.25-field1.4E9-CEP-0.5-theta0-SinFourth-SinPower2-2/lambda745.61e-9-tau33.55e-15-field7E7-CEP0.0-delay0e-15-SinFourth/td.general/laser"
CEPm03="./lambda2237e-9-tau33.55e-15-dt0.25-field1.4E9-CEP-0.3125-theta0-SinFourth-SinPower2-2/lambda745.61e-9-tau33.55e-15-field7E7-CEP0.0-delay0e-15-SinFourth/td.general/laser"
CEPm025="./lambda2237e-9-tau33.55e-15-dt0.25-field1.4E9-CEP-0.25-theta0-SinFourth-SinPower2-2/lambda745.61e-9-tau33.55e-15-field7E7-CEP0.0-delay0e-15-SinFourth/td.general/laser"
CEPm0375="./lambda2237e-9-tau33.55e-15-dt0.25-field1.4E9-CEP-0.375-theta0-SinFourth-SinPower2-2/lambda745.61e-9-tau33.55e-15-field7E7-CEP0.0-delay0e-15-SinFourth/td.general/laser"
CEP025="./lambda2237e-9-tau33.55e-15-dt0.25-field1.4E9-CEP0.25-theta0-SinFourth-SinPower2-2/lambda745.61e-9-tau33.55e-15-field7E7-CEP0.0-delay0e-15-SinFourth/td.general/laser"
CEPm043="./lambda2237e-9-tau33.55e-15-dt0.25-field1.4E9-CEP-0.4375-theta0-SinFourth-SinPower2-2/lambda745.61e-9-tau33.55e-15-field7E7-CEP0.0-delay0e-15-SinFourth/td.general/laser"
CEP05="./lambda2237e-9-tau33.55e-15-dt0.25-field1.4E9-CEP0.5-theta0-SinFourth-SinPower2-2/lambda745.61e-9-tau33.55e-15-field7E7-CEP0.0-delay0e-15-SinFourth/td.general/laser"

# p CEPm05 u 2:3 w l t '-0.5', \

# p CEPm043 u 2:3 w l t '-0.43', \
p CEPm0375 u 2:3 w l t '-0.375'
# CEPm03 u 2:3 w l t '-0.3'
# CEPm025 u 2:3 w l t '-0.25'
# CEP0 u 2:3 w l t '0', \
# CEP025 u 2:3 w l t '0.25'
# CEP05 u 2:3 w l t '0.5"
