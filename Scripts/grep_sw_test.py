##### A script to test the parser scripts, the same reads as parsed by the script
# and by hand are compared, if they are the same and some statistics are returned



# file_by_hand = "/home/luisa/Documents/Work/Uni_Cambridge/Data/Motifs_by_hand_plasmid_mix1_rep1_dil128in1000"
# file_by_script = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/script_qc_aa"
#
# file_by_hand = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/20_singletons/Motifs_by_hand_20_singleton_seqs.fasta"
# file_by_script = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/20_singletons/Motifs_by_script_20_singleton_seq_20_singletons_aa"
#
# file_by_hand = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/20_singletons/Motifs_by_hand_20_singletons_no_insert"
# file_by_script = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/20_singletons/aa_20_singletons_no_inserts"
#
#
#
# file_by_hand = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/quality_control/Motifs_by_hand_5_singletons"
# file_by_script = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/quality_control/aa_5_singletons"
#
# file_by_hand = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/20_singletons/Motifs_by_hand_20_singletons_no_insert"
# file_by_script = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/20_singletons_output/20_singletons_no_insert_aa"

file_by_hand = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/semiglobal_test_run/semiglobal20uniques/20_uniques_by_hand"
file_by_script = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/semiglobal_test_run/semiglobal20uniques/20_uniques_aa"

by_hand = {}
with open(file_by_hand, 'r') as file:
    for i, line in enumerate(file):
        if i%3 == 0:
            seq_id = line.strip()
        if i%3 == 2:
            seq = line.strip()
            by_hand[seq_id] = seq


by_script = {}
with open(file_by_script, 'r') as file:
    for i, line in enumerate(file):
        if i%2 == 0:
            seq_id = line.strip()
        if i%2 == 1:
            seq = ''.join((line.strip()).split(' '))
            by_script[seq_id] = seq




# is there in both and script == hand
# is there in both and script != hand
# is not there in script but by hand
# is not there in both
# is there in script but not by hand


in_both_and_equal = 0
in_both_not_equal = 0
in_script_not_hand = 0
in_hand_not_script = 0
in_neither = 0
total = len(by_hand.keys())

for id in by_hand.keys():
    if id in by_script.keys():
        if by_hand[id] == by_script[id]:
            in_both_and_equal +=1
        elif by_hand[id] == '-':
            in_script_not_hand +=1
            print(f"Script found, Hand not: {id}")
        else:
            in_both_not_equal +=1
            print(f"Different in Hand and script {id}")
    elif by_hand[id] == '-':
        in_neither += 1
    else:
        in_hand_not_script += 1
        print(f"found in hand but not script {id}")
        # print(by_hand[id])

print("\n\n\n\n")
print(f"% sequences found by both which are the same\t\t{in_both_and_equal/total * 100}%")
print(f"% sequences found by both which are not the same\t{in_both_not_equal/total * 100}%")
print(f"% sequences found by hand and not the script\t\t{in_hand_not_script/total * 100}%")
print(f"% sequences found by script and not by hands\t\t{in_script_not_hand/total * 100}%")
print(f"% sequences found by neither\t\t\t\t\t\t{in_neither/total * 100}%")

