from collections import OrderedDict

alphabets = 'abcdefghijklmnopqrstuvwxyz'.upper()
indexes = [i for i in range(len(alphabets))]
positioning_dict = OrderedDict(zip(alphabets, indexes))


def letter_index(lettered_group, positioning_dict):
    index = positioning_dict[lettered_group[0]]
    return index


def reorder_groups(groups):
    letter_order = []
    for i in groups:
        index = positioning_dict[i[0]]
        letter_order.append(index)
    trial = list(zip(groups, letter_order))
    trial.sort(key=lambda tup: tup[1])
    reordered_groups = [i[0] for i in trial]

    return reordered_groups


def Simplify(group_list, positioning_dict, alphabets):
    groups = []
    indexes = []
    for i in group_list:
        index = letter_index(i, positioning_dict)
        indexes.append(index)
    anchor = indexes[0]

    new_indexes = [i - anchor for i in indexes]
    new_letters = [alphabets[i] for i in new_indexes]
    naked = [i[1:] for i in group_list]

    compact = [new_letters[i] + naked[i] for i in range(len(new_letters))]
    compact[0] = 'R' + compact[0][1:]
    compact = " ".join(compact)
    return compact


def Iterator(og):
    split_list = og.split(' ')
    #reordered = reorder_groups(split_list)
    reordered = split_list

    if not reordered:
        new_wln = 'RH'
        return new_wln

    for a in range(len(reordered)):
        if len(reordered[a]) == 1:
            reordered[a] = reordered[a] + '1'

    if len(reordered) > 1:
        new_wln = Simplify(reordered, positioning_dict, alphabets)
        return new_wln

    else:
        new_wln = 'R' + reordered[0][1:]
        return new_wln