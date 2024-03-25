from rest_simulation import REST_simulation

# Class initialization and definitions
apo = REST_simulation(name='apo',
                      rest_folder='/home/s.emanuele/AR/R2_R3_69aa_MD/start_REST/',
                      sequence_offset=378,
                      number_of_replicas=20,
                      output_folder='/home/s.emanuele/REST-Analysis/new_AR_apo/')


l_67q = REST_simulation(name='NT1-67Q',
                         rest_folder='/home/s.emanuele/AR/R2_R3_69aa_NT1Q_MD/drug_confs/start_REST/',
                         sequence_offset=378,
                         number_of_replicas=20,
                         output_folder='/home/s.emanuele/REST-Analysis/new_AR_NT1-67Q/')

l_67q_REST3_05 = REST_simulation(name='NT1-67Q_REST3_05',
                         rest_folder='/home/s.emanuele/AR/R2_R3_69aa_NT1-67Q_REST3_05/REST/',
                         sequence_offset=378,
                         number_of_replicas=20,
                         output_folder='/home/s.emanuele/REST-Analysis/new_AR_NT1-67Q_REST3_05/')

l_69nc = REST_simulation(name='NT1-69NC',
                         rest_folder='/home/s.emanuele/AR/R2_R3_69aa_NT1NC_MD_new/drug_confs/REST/',
                         sequence_offset=378,
                         number_of_replicas=20,
                         output_folder='/home/s.emanuele/REST-Analysis/new_AR_NT1-69NC/')

l_262 = REST_simulation(name='NT1-262',
                        rest_folder='/home/s.emanuele/AR/R2_R3_69aa_NT1-262/REST/',
                        sequence_offset=378,
                        number_of_replicas=20,
                        output_folder='/home/s.emanuele/REST-Analysis/new_AR_NT1-262/')
                    
l_333 = REST_simulation(name='NT1-333',
                        rest_folder='/home/s.emanuele/AR/R2_R3_69aa_NT1-333/REST/',
                        sequence_offset=378,
                        number_of_replicas=20,
                        output_folder='/home/s.emanuele/REST-Analysis/new_AR_NT1-333/')

l_249 = REST_simulation(name='NT1-249',
                        rest_folder='/home/s.emanuele/AR/R2_R3_69aa_NT1-249/REST/',
                        sequence_offset=378, # The first H is 380. But we also start with a CAP, so 378 instead of 379.
                        number_of_replicas=20,
                        output_folder='/home/s.emanuele/REST-Analysis/new_AR_NT1-249/')
                    
l_269 = REST_simulation(name='NT1-269',
                        rest_folder='/home/s.emanuele/AR/R2_R3_69aa_NT1-269/REST/',
                        sequence_offset=378, # The first H is 380. But we also start with a CAP, so 378 instead of 379.
                        number_of_replicas=20,
                        output_folder='/home/s.emanuele/REST-Analysis/new_AR_NT1-269/')

l_maso = REST_simulation(name='Masofaniten',
                        rest_folder='/home/s.emanuele/AR/R2_R3_69aa_Masofaniten/REST_Jaqi/',
                        sequence_offset=389, # The first H is 380. But we also start with a CAP, so 378 instead of 379.
                        number_of_replicas=16,
                        output_folder='/home/s.emanuele/REST-Analysis/Masofaniten/')

l_1aa = REST_simulation(name='1AA',
                        rest_folder='/home/s.emanuele/AR/R2_R3_69aa_1AA/REST_Jiaqi/',
                        sequence_offset=389, # The first H is 380. But we also start with a CAP, so 378 instead of 379.
                        number_of_replicas=16,
                        output_folder='/home/s.emanuele/REST-Analysis/1AA/')


to_analyse = [apo, l_67q, l_67q_REST3_05, l_69nc, l_333, l_249, l_269, l_maso, l_1aa]
# to_streamlit = [apo, l_67q, l_69nc, l_262, l_333, l_249, l_269, l_maso, l_1aa]
to_streamlit = [l_maso, l_1aa]
