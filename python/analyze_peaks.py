#! /usr/bin/python

import math
import sys

def parse_peak_data(f):
  res = []
  for line in f:
    line.strip()
    parts = line.split(' ')
    res.append((int(parts[0]), int(parts[1]), float(parts[2])))
  return res


def compute_l2_distances(peaks1, peaks2):
  l2_dsts = []
  for p1, p2 in zip(peaks1, peaks2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    d = math.sqrt(dx * dx + dy * dy)
    l2_dsts.append(d)
  return l2_dsts


def mean(l):
  s = 0.0
  for x in l:
    s += x
  return s / len(l)


def median(l):
  l2 = sorted(l)
  return l2[len(l) / 2]


if __name__ == '__main__':
  name1 = sys.argv[1]
  name2 = sys.argv[2]
  with open(name1, 'r') as f1:
    with open(name2, 'r') as f2:
      peaks1 = parse_peak_data(f1)
      peaks2 = parse_peak_data(f2)
      l2_dsts = compute_l2_distances(peaks1, peaks2)
      print 'min: {}  max: {}  mean: {}  median: {}'.format(min(l2_dsts),
                                                            max(l2_dsts),
                                                            mean(l2_dsts),
                                                            median(l2_dsts))
