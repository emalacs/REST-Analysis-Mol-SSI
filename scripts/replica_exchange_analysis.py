import os
from util import run_script, run_script_sp, file_exists, create_folder, set_working_directory
import MDAnalysis as mda
from MDAnalysis import transformations
import pandas as pd


def make_whole(replica_directory, contains_ligand, gmx_gpu_mpi_path, prefix):
    print(f'Making whole {replica_directory}')

    if contains_ligand:
        gmx_convert_cmd = f'echo 18 | {gmx_gpu_mpi_path} convert-tpr -s {replica_directory}/{prefix}.tpr -o {replica_directory}/protein.tpr'
    else:
        gmx_convert_cmd = f'echo 14 | {gmx_gpu_mpi_path} convert-tpr -s {replica_directory}/{prefix}.tpr -o {replica_directory}/protein.tpr'
    run_script_sp(gmx_convert_cmd)
    file_exists(f'{replica_directory}/protein.tpr')
    run_script_sp(f'{gmx_gpu_mpi_path} editconf -f {replica_directory}/protein.tpr -o {replica_directory}/protein.gro')
    file_exists(f'{replica_directory}/protein.gro')
    run_script_sp(f'echo 0 | {gmx_gpu_mpi_path} trjconv -s {replica_directory}/protein.tpr -f {replica_directory}/{prefix}.xtc -o {replica_directory}/whole.xtc -pbc whole')
    file_exists(f'{replica_directory}/whole.xtc')
    print(f'Replica {replica_directory} is whole again')


def center_wrap(replica, format='xtc'):
    '''
    This function is the same from pbc4solutes.py
    '''

    print(f'Centering replica {replica}')
    topology = f'{replica}/protein.tpr'
    trajectory = f'{replica}/whole.xtc'

    u = mda.Universe(topology, trajectory)
    system = u.select_atoms('all')
    protein = u.select_atoms('protein')
    not_protein = u.select_atoms('not protein')
    transforms = [transformations.unwrap(system), transformations.center_in_box(protein, wrap=False), transformations.wrap(not_protein, compound='residues')]
    u.trajectory.add_transformations(*transforms)
    with mda.Writer(f'{replica}/pbc_corrected.'+format, system.n_atoms) as W:
        for ts in u.trajectory:
            W.write(system)
    print(f'Replica {replica} centered.')


def make_demux(simulation_directory, n_replicas, gmx_gpu_mpi_path, prefix, simulation_time_step = 0.002):
    '''
    By default the demux reading is done by reading the 0 replica log.
    '''
    
    working_directory = os.getcwd()
    demux_directory = f'{simulation_directory}/demux/'
    create_folder(demux_directory)
    create_folder(f'{demux_directory}/repex_analysis/')

    set_working_directory(demux_directory)
    with open('demux.fix.pl', 'w') as file:
        file.write(demux_fix_pl)

    if not os.path.exists(f'{demux_directory}/demux.fix.pl'):
        print(f'{demux_directory}/demux.fix.pl does not exist, check the folder.')
        exit()
    # if not os.path.exists(f'{demux_directory}/make_demux.sh'):
    #     print('make_demux.sh does not exist, check the folder.')
    #     exit()
    # if not os.path.exists(f'{demux_directory}/repex_analysis.sh'):
    #     print('repex_analysis.sh does not exist, check the folder.')
    #     exit()
    

    print(f'Set demux working directory as {demux_directory}.\nStarting demux analysis')#.\nAll files will be saved in {out_demux_directory}')
    # demux_cmd = f'echo {simulation_time_step} | perl demux.fix.pl ../{n_replicas[0]}/production.log'
    demux_cmd = f'echo {simulation_time_step} | perl demux.fix.pl ../0/{prefix}.log'
    print(f'Running demux:\n\t {demux_cmd}')
    run_script_sp(demux_cmd)
    if not os.path.exists(f'{demux_directory}/replica_index.xvg'):
        print(f'replica_index.xvg does not exist. Something went wrong on {demux_cmd}.')
        exit()
    if not os.path.exists(f'{demux_directory}/replica_temp.xvg'):
        print(f'replica_temp.xvg does not exist. Something went wrong on {demux_cmd}.')
        exit()


    # get_demux_trajectories(simulation_directory, demux_directory, n_replicas, out_demux_directory)
    get_demux_trajectories(simulation_directory, demux_directory, n_replicas, gmx_gpu_mpi_path)
    # run_script(f'cp {demux_directory}/replica_index.xvg {out_demux_directory}')
    # run_script(f'cp {demux_directory}/replica_temp.xvg {out_demux_directory}')
    print(f'demux analysis completed, check files')
    set_working_directory(working_directory)


# def get_demux_trajectories(simulation_directory, demux_directory, n_replicas, out_demux_directory):
def get_demux_trajectories(simulation_directory, demux_directory, n_replicas, gmx_gpu_mpi_path):
    reps = f'{n_replicas[0]}..{n_replicas[-1]}'
    #awk_cmd_trr = f"awk '{{if($1==0){{print}} if(n==500){{$1=$1-800.0; print;n=0}} n++;}}' {demux_directory}/replica_index.xvg > {demux_directory}/replica_index.n500.s0.-800.xvg"
    #trjcat_cmd_trr = f'gmx trjcat -f {simulation_directory}/{{{reps}}}/production.trr -demux {demux_directory}/replica_index.n500.s0.-800.xvg -o {demux_directory}/{{{reps}}}.trr'
    #run_script(awk_cmd_trr)
    #run_script(trjcat_cmd_trr)
    # if not os.path.exists(f'{demux_directory}/replica_index.n500.s0.-800.xvg'):
            # print(f'replica_index.n500.s0.-800.xvg does not exist. Something went wront on {awk_cmd_trr}.')
            # exit()
    #run_script(f'cp {demux_directory}/replica_index.n500.s0.-800.xvg {out_demux_directory}')
    awk_cmd_xtc = f"awk '{{if($1==0){{print}} if(n==50){{$1=$1-80.0; print;n=0}} n++;}}' {demux_directory}/replica_index.xvg > {demux_directory}/replica_index.n50.s0.-80.xvg"
    run_script_sp(awk_cmd_xtc)
    trjcat_cmd_xtc = f'{gmx_gpu_mpi_path} trjcat -f {simulation_directory}/{{{reps}}}/pbc_corrected.xtc -demux {demux_directory}/replica_index.n50.s0.-80.xvg -o {demux_directory}/{{{reps}}}.xtc'
    run_script_sp(trjcat_cmd_xtc)
    # run_script(f'cp {demux_directory}/replica_index.n50.s0.-80.xvg {out_demux_directory}')

    if not os.path.exists(f'{demux_directory}/replica_index.n50.s0.-80.xvg'):
            print(f'replica_index.n500.s0.-800.xvg does not exist. Something went wront on {awk_cmd_xtc}.')
            exit()


def get_acceptance_ratios(rest_folder, prefix):

    working_directory = os.getcwd()
    set_working_directory(rest_folder)
    with open('AcceptRatio20.sh', 'w') as file:
        file.write(write_acceptratio20_sh(prefix))
    run_script_sp('bash AcceptRatio20.sh > acceptance_ratios.txt')
    # run_script(f'cp acceptance_ratios.txt {out_demux_directory}')
    acceptance_ratios_df = pd.read_csv('acceptance_ratios.txt', header=None)
    acceptance_ratios_df = acceptance_ratios_df[~acceptance_ratios_df[0].str.startswith('Repl pr')]
    acceptance_ratios_df.columns = ['values']
    mask = acceptance_ratios_df['values'].str.contains('x')
    temp_exchange_df = acceptance_ratios_df[mask].copy()
    temp_exchange_df.columns = ['exchanging_replicas']
    temp_exchange_df.reset_index(inplace=True)
    temp_exchange_df.drop(columns=['index'], inplace=True)
    temp_values_df = acceptance_ratios_df[~mask].copy()
    temp_values_df = temp_values_df.iloc[1:]
    temp_values_df[['succesfull_exchanges', 'total', 'acceptance_ratio']] = temp_values_df['values'].str.split(' ', expand=True)
    temp_values_df.drop(columns=['values'], inplace=True)
    temp_values_df.reset_index(inplace=True)
    temp_values_df.drop(columns=['index'], inplace=True)

    acceptance_ratios_df = pd.concat([temp_exchange_df, temp_values_df], axis=1)
    acceptance_ratios_df.columns = ['Exchanging Replicas', 'Successful Exchanges', 'Total', 'Acceptance Ratio']

    del temp_values_df, temp_exchange_df
    set_working_directory(working_directory)
    
    return acceptance_ratios_df


def get_repex_analysis(rest_folder, n_replicas, prefix):

    working_directory = os.getcwd()
    set_working_directory(f'{rest_folder}/demux/')
    with open('repex_analysis.sh', 'w') as file:
        file.write(write_repex_analysis_sh(prefix))
    run_script_sp(f'bash repex_analysis.sh')

    colnames = ['time'] + n_replicas
    replica_temperature_df = pd.read_csv(f'{rest_folder}/demux/replica_temp.xvg', header=None, sep='\s+', names=colnames)
    replica_temperature_df['demux'] = 'temperature'
    replica_index_df = pd.read_csv(f'{rest_folder}/demux/replica_index.n50.s0.-80.xvg', header=None, sep='\s+', names=colnames)
    replica_index_df['demux'] = 'demux'
    repex_analysis_df = pd.concat([replica_temperature_df, replica_index_df])


    # TODO check if those files are good as well compared to the ones I am using now, from the original code.
    # repex_analysis_df = pd.DataFrame()
    # for replica in n_replicas:
    #     temp_df = pd.read_csv(f'repex_analysis/rep.temp.{replica}.dat', header=None, )
    #     temp_df.columns=[f'temperature_{replica}']
    #     repex_analysis_df = pd.concat([repex_analysis_df, temp_df], axis=1)    
    #     temp_df = pd.read_csv(f'repex_analysis/rep.index.{replica}.dat', header=None, )
    #     temp_df.columns=[f'index_{replica}']
    #     repex_analysis_df = pd.concat([repex_analysis_df, temp_df], axis=1)

    set_working_directory(working_directory)

    return repex_analysis_df


def get_round_trip(rest_folder):
    
    working_directory = os.getcwd()
    set_working_directory(f'{rest_folder}/demux/')
    with open('round_trip.sh', 'w') as file:
        file.write(round_trip_sh)
    run_script_sp(f'bash round_trip.sh')
    round_trip_df = pd.read_csv('RoundTrip.stats.dat', header=None)
    round_trip_df.columns = ['values']

    round_trip_df[['Replica', 'Round Trip']] = round_trip_df['values'].str.split(' ', expand=True)
    round_trip_df.set_index('Replica', inplace=True)
    round_trip_df.drop(columns=['values'], inplace=True)
    set_working_directory(working_directory)

    return round_trip_df

# SCRIPTS
demux_fix_pl = '''
#!/usr/bin/perl -w

# in: input filename
$in = shift || die("Please specify input filename");
# If your exchange was every N ps and you saved every M ps you can make for
# the missing frames by setting extra to (N/M - 1). If N/M is not integer,
# you're out of luck and you will not be able to demux your trajectories at all.
$extra = shift || 0 ;
$ndx  = "replica_index.xvg";
$temp = "replica_temp.xvg";

@comm = ("-----------------------------------------------------------------",
	 "Going to read a file containing the exchange information from",
	 "your mdrun log file ($in).", 
	 "This will produce a file ($ndx) suitable for",
	 "demultiplexing your trajectories using trjcat,",
	 "as well as a replica temperature file ($temp).",
	 "Each entry in the log file will be copied $extra times.",
	 "-----------------------------------------------------------------");
for($c=0; ($c<=$#comm); $c++) {
    printf("$comm[$c]\\n");
}

print "Time step of the simulation (ps)? ";
$dtime=<>;

# Open input and output files
open (IN_FILE,"$in") || die ("Cannot open input file $in");
open (NDX,">$ndx") || die("Opening $ndx for writing");
open (TEMP,">$temp") || die("Opening $temp for writing");


sub pr_order {
    my $t     = shift;
    my $nrepl = shift;
    printf(NDX "%-20.2f",$t);
    for(my $k=0; ($k<$nrepl); $k++) {
	my $oo = shift;
	printf(NDX "  %3d",$oo);
    }
    printf(NDX "\\n");
}

sub pr_revorder {
    my $t     = shift;
    my $nrepl = shift;
    printf(TEMP "%-20.2f",$t);
    for(my $k=0; ($k<$nrepl); $k++) {
	my $oo = shift;
	printf(TEMP "  %3d",$oo);
    }
    printf(TEMP "\\n");
}

$nrepl = 0;
$init  = 0;
$tstep = 0;
$nline = 0;
$tinit = 0;
while ($line = <IN_FILE>) {
    chomp($line);
    
    if (index($line,"init_t") >= 0) {
	@log_line = split (' ',$line);
	$tinit = $log_line[2];
    }
    if (index($line,"Repl") == 0) {
	@log_line = split (' ',$line);
	if (index($line,"There") >= 0) {
	    $nrepl = $log_line[3];
	}
	elsif (index($line,"time") >= 0) {
	    $tstep = ($dtime)*$log_line[4];
	}
	elsif ((index($line,"Repl ex") == 0) && ($nrepl == 0)) {
            # Determine number of replicas from the exchange information
	    printf("%s\\n%s\\n",
		   "WARNING: I did not find a statement about number of replicas",
		   "I will try to determine it from the exchange information.");
	    for($k=2; ($k<=$#log_line); $k++) {
		if ($log_line[$k] ne "x") {
		    $nrepl++;
		}
	    }
	}
	if (($init == 0) && ($nrepl > 0)) {
	    printf("There are $nrepl replicas.\\n");

	    @order = ();
            @revorder = ();
	    for($k=0; ($k<$nrepl); $k++) {
		$order[$k] = $k;
                $revorder[$k] = $k;
	    }
	    for($ee=0; ($ee<=$extra); $ee++) {
		pr_order($tinit+$ee,$nrepl,@order);
		pr_revorder($tinit+$ee,$nrepl,@revorder);
		$nline++;
	    }
	    $init = 1;
	}

	if (index($line,"Repl ex") == 0) {
	    $k = 0;
	    for($m=3; ($m<$#log_line); $m++) {
		if ($log_line[$m] eq "x") {
                    $mp1 = $log_line[$m+1];
                    $mm1 = $log_line[$m-1];
		    $revorder[$order[$mm1]] = $mp1;
		    $revorder[$order[$mp1]] = $mm1;
		    $tmp = $order[$mm1];
		    $order[$mm1] = $order[$mp1];
		    $order[$mp1] = $tmp;
#	    printf ("Swapping %d and %d on line %d\\n",$k,$k+1,$line_number); 
		}
		else {
		    $k++;
		}
	    }
	    for($ee=0; ($ee<=$extra); $ee++) {
		pr_order($tstep+$ee,$nrepl,@order);
		pr_revorder($tstep+$ee,$nrepl,@revorder);
		$nline++;
	    }
	}
    }
}
close IN_FILE;
close NDX;
close TEMP;

printf ("Finished writing $ndx and $temp with %d lines\\n",$nline);

'''

def write_acceptratio20_sh(prefix):
    acceptratio20_sh = f'''
grep "Repl pr" 0/{prefix}.log
#off by factor of 2 - not splitting exchanges
n=`grep "Repl ex" 0/{prefix}.log | wc -l`
echo $n
for i in {{0..8}}
do 
j=$[$i+1]
echo "$i x  $j"
a=`grep "Repl ex" $i/{prefix}.log | grep "$i x  $j" | wc -l`
echo $a $n > ex
awk '{{print($1,$2/2,$1/$2*2)}}' ex
done

for i in {{9..14}}
do
j=$[$i+1]
echo "$i x $j"
#grep "Repl ex" $i/{prefix}.log | wc -l
#grep "Repl ex" $i/{prefix}.log | grep "$i x  $j" | wc -l 
a=`grep "Repl ex" $i/{prefix}.log | grep "$i x $j" | wc -l`
echo $a $n > ex
awk '{{print($1,$2/2,$1/$2*2)}}' ex
done

for i in {{15..18}}
do
    j=$[$i+1]
    echo "$i x $j"
    #grep "Repl ex" $i/{prefix}.log | wc -l
    #grep "Repl ex" $i/{prefix}.log | grep "$i x  $j" | wc -l 
    a=`grep "Repl ex" $i/{prefix}.log | grep "$i x $j" | wc -l`
    echo $a $n > ex
    awk '{{print($1,$2/2,$1/$2*2)}}' ex
done

'''
    return acceptratio20_sh


def write_repex_analysis_sh(prefix):
    repex_analysis_sh = f'''

labhome=/dartfs-hpc/rc/lab/R/RobustelliP/
#demux.pl 1/{prefix}.log > demux.out
#mkdir repex_anal
#reps=`grep "There are " demux.out | awk '{{print($3)}}'` 
#k=$(($reps))
#echo $k
# mkdir repex_analysis
for i in {{0..19}}
do
echo $i
j=$[$i+2]
awk '{{print($'$j')}}' replica_temp.xvg > repex_analysis/rep.temp.$i.dat
awk '{{print($'$j')}}' replica_index.xvg > repex_analysis/rep.index.$i.dat
# nohup python3 1d.timeseries.hist.py repex_analysis/rep.temp.$i.dat repex_analysis/rep.temp.$i &
# nohup python3 1d.timeseries.hist.plotly.py repex_analysis/rep.temp.$i.dat repex_analysis/rep.temp.$i.plotly $i &
# wait
# nohup python3 1d.timeseries.hist.py repex_analysis/rep.index.$i.dat repex_analysis/rep.index.$i &
# nohup python3 1d.timeseries.hist.plotly.py repex_analysis/rep.index.$i.dat repex_analysis/rep.index.$i.plotly $i &
# wait

done

'''
    return repex_analysis_sh

round_trip_sh = f'''

>RoundTrip.stats.dat
for i in {{0..19}}
do
awk '{{n++; if($1==19){{d=1}} if($1==0) if(d==1){{d=0; rt+=1;}} if(rt>0) print("'$i'",n/rt*400*.002/1000)}}' repex_analysis/rep.temp.$i.dat | tail -n1 >> RoundTrip.stats.dat
done

awk '{{a+=$2; n++; print("Average:",a/n)}}' RoundTrip.stats.dat | tail -n1 > RT.ave
paste RT.ave >> RoundTrip.stats.dat
paste RoundTrip.stats.dat

'''