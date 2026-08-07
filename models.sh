# sweep a no-drug, no-plasmid community and output the resulting trajectories
octave-cli model/octave_compat/run_sweep.m model/params/params.par --no-drug --tmax 1100 --dt 4 --solver ode45 --no-mse --no-summary --sweep init_As=log:1e-4:10:6 --sweep init_Bs=log:1e-4:10:6 --sweep init_Cs=log:1e-4:10:6 --traj data/model/inoculum_sweep.tsv
# single run with no drug and no plasmid, same initial inoculum as 1E8
octave-cli model/octave_compat/run_sim.m model/params/params.par --no-drug --tmax 1100 --init As=10 --init Bs=10 --init Cs=10 --out data/model/nodrug.tsv
# sweep drug efficacy parameters in a no-plasmid community
octave-cli model/octave_compat/run_sweep.m model/params/params.par --drug --no-conjugation --tmax 1100 --dt 4 --solver lsode --no-mse --no-summary --init As=10 --init Bs=10 --init Cs=10 --sweep Ab_init=0,0.0001,0.000215,0.000464,0.001,0.00215,0.00464,0.01,0.0215,0.0464,0.1 --sweep gamma=log:0.003:3:16 --traj data/model/ab_gamma_no_plasmid.tsv
# sweep epsilon alone in a no-plasmid community, at the same dose as the runs below
octave-cli model/octave_compat/run_sweep.m model/params/params.par --drug --no-conjugation --tmax 1100 --dt 4 --solver lsode --no-mse --no-summary --ab0 0.01 --gamma 0.019 --init As=10 --init Bs=10 --init Cs=10 --sweep epsilon=log:0.01:100:21 --traj data/model/epsilon_all_sensitive.tsv --traj-vars As,Bs,Cs
# sweep epsilon and eta in a LM-plasmid community
octave-cli model/octave_compat/run_sweep.m model/params/params.par --community 3 --tmax 1100 --dt 4 --solver lsode --no-mse --no-summary --ab0 0.01 --gamma 0.019 --init As=0 --init Ar=10 --init Bs=10 --init Cs=10 --sweep epsilon=log:0.01:100:21 --sweep eta=log:1e-4:1e-2:16 --traj data/model/epsilon_eta_plasmid.tsv --traj-vars Ar,Bs,Br,Cs,Cr
# sweep epsilon in the non-transmissible community
octave-cli model/octave_compat/run_sweep.m model/params/params.par --community 4 --tmax 1100 --dt 4 --solver lsode --no-mse --no-summary --ab0 0.01 --gamma 0.019 --init As=0 --init Ar=10 --init Bs=10 --init Cs=10 --sweep epsilon=log:0.01:100:21 --traj data/model/epsilon_no_transfer.tsv --traj-vars Ar,Bs,Br,Cs,Cr
# rescue: conjugation just outruns killing, PM/PL are converted instead of lost.
octave-cli model/octave_compat/run_sim.m model/params/params.par --community 3 --tmax 600 --dt 2 --solver lsode --ab0 0.01 --gamma 0.019 --epsilon 1.584893192 --eta 0.0008576958986 --init As=0 --init Ar=10 --init Bs=10 --init Cs=10 --out data/model/exp_match_PVI.tsv
# same drug, plasmid cannot transfer (community 4 forces eta = 0): washes out
octave-cli model/octave_compat/run_sim.m model/params/params.par --community 4 --tmax 600 --dt 2 --solver lsode --ab0 0.01 --gamma 0.019 --epsilon 1.584893192 --init As=0 --init Ar=10 --init Bs=10 --init Cs=10 --out data/model/exp_match_PVB.tsv
# conjugation works but is too slow
octave-cli model/octave_compat/run_sim.m model/params/params.par --community 3 --tmax 600 --dt 2 --solver lsode --ab0 0.01 --gamma 0.019 --epsilon 1.584893192 --eta 0.0004641588834 --init As=0 --init Ar=10 --init Bs=10 --init Cs=10 --out data/model/exp_match_slow_conj.tsv
# same conjugation as the rescue, one epsilon step up: killing outruns transfer
octave-cli model/octave_compat/run_sim.m model/params/params.par --community 3 --tmax 600 --dt 2 --solver lsode --ab0 0.01 --gamma 0.019 --epsilon 2.511886432 --eta 0.0008576958986 --init As=0 --init Ar=10 --init Bs=10 --init Cs=10 --out data/model/exp_match_strong_drug.tsv
# no plasmid anywhere: same drug and same conjugation rate, but LM starts sensitive
octave-cli model/octave_compat/run_sim.m model/params/params.par --community 3 --tmax 600 --dt 2 --solver lsode --ab0 0.01 --gamma 0.019 --epsilon 1.584893192 --eta 0.0008576958986 --init As=10 --init Ar=0 --init Bs=10 --init Cs=10 --out data/model/exp_match_no_plasmid.tsv
