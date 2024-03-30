import time

return_mols_per_file = []
for j in range(3):
    ct = 0
    for i in range(len(return_mols_per_file)):
        ct += len(return_mols_per_file[i][j])
    print(ct)

for k in range(3):
    t = time.time()
    for i in range(len(return_mols_per_file)):
        for j in range(i + 1, len(return_mols_per_file)):
            for keys in return_mols_per_file[i][k].keys():
                if keys in return_mols_per_file[j][k]:
                    return_mols_per_file[j][k].pop(keys)
    print(time.time() - t)

for j in range(3):
    ct = 0
    for i in range(len(return_mols_per_file)):
        ct += len(return_mols_per_file[i][j])
    print(ct)

train = {}
valid = {}
test = {}
for j in range(3):
    for i in range(len(return_mols_per_file)):
        for keys in return_mols_per_file[i][j]:
            if j == 0:
                train[keys] = 0
            elif j == 1:
                valid[keys] = 0
            elif j == 2:
                test[keys] = 0
