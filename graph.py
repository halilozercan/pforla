import pickle
import random
import sys

import editdistance as editdistance
from gatb import Bank, Graph

from pathgraph import PathGraph

halt_prob = 0.02


def solid_kmers(sequence, graph, kmer_size):
    result = []
    sit = set()
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence.sequence[i:(i + kmer_size)]
        node = graph[kmer]
        if node in graph and node.out_degree > 0 and node.in_degree > 0:
            result.append((node, i))
            sit.add(i)
            # print(node)
            # print(node.succs)

    """
    minhash = MinHash([str(i) for i in nodes], 4, sim=0.7)

    for i in range(len(sequence) - kmer_size + 1):
        if i in sit:
            continue
        kmer = sequence.sequence[i:(i + kmer_size)]
        results = minhash.query(kmer)
        for candidate_kmer in results:
            node = graph[candidate_kmer]
            if node.out_degree > 0:
                result.append((node, i))
    """
    result = list(set(result))
    return sorted(result, key=lambda x: x[1])


def find_path(seed, target, ref_seq, min_editdistance=10000000):
    stack = [(seed, (i for i in seed.succs), str(seed))]
    # editdist = Distance('', '')
    branch = 200
    best_path = seed
    path_found = False
    while stack and branch:
        (vertex, g, path) = stack[-1]
        if random.random() < halt_prob:
            stack.pop()
            continue
        try:
            n = next(g)
        except StopIteration:
            stack.pop()
            branch -= 1
            continue

        path += str(n)[-1]

        # editdist.update(ref_seq[:len(path)], path)
        # distance = editdist.eval()
        distance = editdistance.eval(ref_seq[:len(path)], path)

        if distance < len(path) * 0.4 and len(path) <= len(ref_seq):
            if len(best_path) < len(path) or distance < min_editdistance:
                min_editdistance = distance
                best_path = path
        elif distance > len(path) * 0.4:
            branch -= 1
            stack.pop()
            continue

        if str(n) == str(target):
            path_found = True
            branch -= 1
            #            print("Reporting a path")
            yield path, distance, 0
        else:
            stack.append((n, (i for i in n.succs), path))
    if min_editdistance < 10000 and not path_found:
        yield best_path, min_editdistance, 1


def path_graph(seq, clusters):
    pathgraph = []
    for cluster_index in range(len(clusters) - 1):
        cluster = clusters[cluster_index]
        for seed_index, seed in enumerate(cluster):
            targets = []
            _cluster_index = cluster_index + 1
            while len(targets) < 5 and _cluster_index < len(clusters):
                targets += clusters[_cluster_index]
                _cluster_index += 1
            targets = targets[:5]

            for target in targets:
                ref_seq = seq.sequence[seed[1]:(target[1] + KMER_SIZE)].decode()
                # print('Trying to go from %s to %s %d %d' % (seed[0], target[0], seed[1], target[1]))
                min_editdistance = 1000000
                min_path = None
                for path, weight, ext in find_path(seed[0], target[0], ref_seq, min_editdistance):
                    if ext and len(path) > KMER_SIZE:
                        # print("Found extension at the length of %d\n%s\n%s" % (len(path), path, str(ref_seq[:len(
                        # path)])))
                        new_target = -1 * (seed[1] + len(path) - KMER_SIZE)
                        edge1 = (seed[1], new_target, path, weight)
                        edge2 = (new_target, target[1], '', abs(target[1] - seed[1] - (len(path) - KMER_SIZE)))
                        pathgraph.append(edge1)
                        pathgraph.append(edge2)
                        # print("Appending to path graph:\n%s\n%s" % (edge1, edge2) )

                    elif weight < min_editdistance:
                        min_editdistance = weight
                        min_path = path

                if min_path is not None:
                    # print('%s %s' % (min_path, ref_seq))
                    # print(min_editdistance)
                    # print(limit)
                    pathgraph.append((seed[1], target[1], min_path, min_editdistance))
                else:
                    pathgraph.append((seed[1], target[1], '', target[1] - seed[1]))
                    pass
    # We reached tail. Now special extension comes in
    cluster = clusters[-1]

    return pathgraph


graph_file = sys.argv[1]
reads_file = sys.argv[2]
out_file = sys.argv[3]
KMER_SIZE = int(sys.argv[4])
first_n_reads = int(sys.argv[5])

bank = Bank(reads_file)
print("File '%s' is of type: %s" % (bank.uri, bank.type))

graph = Graph('-in %s' % graph_file)

nseqs = 0

# pickle.dump(minhash, open('minhash', 'wb'))
# bktree = pickle.load(open('BKTREE_19', 'rb'))
print("Got the nodes")

total_path_found = 0
total_read_length = 0

read_paths = []
for i, seq in enumerate(bank):
    if i > first_n_reads:
        continue
    # 'seq' is of type 'Sequence'.
    # Accessing 'Sequence' internals is done as follows:
    #   sequence header : seq.comment
    #   sequence quality: seq.quality (Fastq only)
    #   sequence letters: seq.sequence
    #   sequence size   : len(seq)
    seqid = seq.comment.decode("utf-8").split(" ")[0]
    print('%d: %s: %d letters' % (i, seqid, len(seq)))
    nseqs += 1

    seeds = solid_kmers(seq, graph, KMER_SIZE)
    print('#solid kmers in sequence: %d' % len(seeds))

    if len(seeds) == 0:
        continue

    clusters = []
    safe = True
    for j in range(len(seeds)):
        if safe:
            clusters.append([])

        clusters[-1] += [seeds[j]]
        last_index = seeds[j][1]

        if j + 1 < len(seeds):
            safe = ((seeds[j + 1][1] - seeds[j][1]) > 1)

    print('Total cluster count %d seed' % len(clusters))
    read_path = path_graph(seq, clusters)

    pathgraph = PathGraph()
    for path in read_path:
        # pathgraph.add_vertex(str(path[0]))
        # pathgraph.add_vertex(str(path[1]))
        pathgraph.add_edge(str(path[0]), str(path[1]), path[3])

    # pathgraph.write_to_file('deneme.out')

    shortest_path = pathgraph.shortest_path(str(seeds[0][1]), str(clusters[-1][:5][-1][1]))
    shortest_path.append(str(seeds[0][1]))
    # print(shortest_path)
    total_weight = 0
    for i in range(len(shortest_path) - 1):
        total_weight += pathgraph.vertices[shortest_path[i]][shortest_path[i + 1]]

    normal_distance = clusters[-1][:5][-1][1] - seeds[0][1]
    total_path_found += (normal_distance - total_weight)

    total_read_length += len(seq)

    print(shortest_path)
    print(total_weight)
    # pprint.pprint(read_path)
    read_paths.append((read_path, shortest_path, total_weight))

print("==============================================================")
print("Total path found %d %d" % (total_path_found, total_read_length))
pickle.dump(read_paths, open(out_file, 'wb'))
