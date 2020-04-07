# EVB
This is an example for EVB calculations in water. It reads in input files and save the EVB profile as a fig.png file and print out the set **alpha** and **H12** parameters together with calculated **free activation erergy** and **solvation energy** 

## Files required :
### one .xvg file for each &lambda; step
For Each &lambda; value, the energy for the initial state(lambda=0, bonded state) and final state(lambda=1, unbonded state) are extacted to a $lambda$.txt file. (Which in turn, were generated by grepping the energies from the GROMACS output.xvg for each &lambda; step)

### fep.xvg
The results from FEP calculation(generated by *gmx bar*)were tidied up(i.e. only the change in energy **DG** in the **Final Results** section) to a fep.xvg file.

### solv.txt 
solvation energy as sum of the solute-solvent electrostatic and LJ interactions. Extracted by *gmx energy*, an average value is used for each lambda step. (It can be and possible to get the solvation energy for each specific step, this will soon be implemented) 

### lambda.xvg
A list of all &lambda; values 
