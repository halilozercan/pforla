import numpy


class Distance:
    def __init__(self, a, b):
        self.a = a
        self.b = b

        self.distances = numpy.zeros(shape=(300, 300))
        for i in range(300):
            self.distances[i, 0] = i
        for i in range(300):
            self.distances[0, i] = i

        for i1, c1 in enumerate(a):
            for i2, c2 in enumerate(b):
                if c1 == c2:
                    self.distances[i1 + 1, i2 + 1] = self.distances[i1, i2]
                else:
                    self.distances[i1 + 1, i2 + 1] = 1 + min(self.distances[i1, i2 + 1], self.distances[i1 + 1, i2],
                                                             self.distances[i1, i2])

    def eval(self):
        return self.distances[len(self.a), len(self.b)]

    def remove(self, len_a, len_b):
        self.a = self.a[:-len_a]
        self.b = self.b[:-len_b]

    def update(self, new_a, new_b):
        if len(new_a) > len(self.a):
            self.concatenate(new_a[len(self.a):], new_b[len(self.b):])
        elif len(new_a) < len(self.a):
            self.remove(len(self.a) - len(new_a), len(self.b) - len(new_b))

    def concatenate(self, a, b):
        last_a = len(self.a)
        last_b = len(self.b)
        self.a += a
        self.b += b

        for i in range(last_a, len(self.a)):
            c1 = self.a[i]
            for j, c2 in enumerate(self.b):
                if c1 == c2:
                    self.distances[i + 1, j + 1] = self.distances[i, j]
                else:
                    self.distances[i + 1, j + 1] = 1 + min(self.distances[i, j + 1], self.distances[i + 1, j],
                                                           self.distances[i, j])

        for j in range(last_b, len(self.b)):
            c2 = self.b[j]
            for i, c1 in enumerate(self.a):
                if c1 == c2:
                    self.distances[i + 1, j + 1] = self.distances[i, j]
                else:
                    self.distances[i + 1, j + 1] = 1 + min(self.distances[i, j + 1], self.distances[i + 1, j],
                                                           self.distances[i, j])


a = Distance('AC', 'AC')
print(a.eval())
a.update('ACTGCTGACAGTTTCACAG', 'ACTGCTGACAGTCACACT')
print(a.eval())
