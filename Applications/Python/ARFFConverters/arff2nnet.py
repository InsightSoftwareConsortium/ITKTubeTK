import csv
import sys
import getopt

if __name__ == '__main__':
  try:
    opts, args = getopt.getopt(sys.argv[1:], 's:')
  except getopt.GetoptError, err:
    print str(err)
    sys.exit(1)

  skiplist = []
  for optname, optarg in opts:
    if optname == '-s':
      try:
        if optarg[0] == '=':
          optarg = optarg[1:]

        skiplist = [int(x) for x in optarg.split(',')]
        def throwerror():
          raise Exception("Problem")
        [throwerror() for x in skiplist if x <0]
      except:
        print 'error converting arguments to skip'
        print optarg
        sys.exit(2)
          
  if len(args) != 1:
    print 'need input'
    sys.exit(1)

  f = open(args[0], 'r')
  for line in f:
    if 'class' in line:
      classes = line[line.find('{')+1:line.find('}')].split(',')
    if '@DATA' in line:
      break

  # start main processin
  for line in f:
    features = line.split(',')
    features = [x for (i, x) in enumerate(features) if not i in skiplist]
    clazz = features[-1].strip()
    features = features[:-1]

    classvector = ['-1'] * 3
    classvector[classes.index(clazz)] = '1'

    features.extend(classvector)

    #features = ['%d:%s' % (i+1,f) for (i,f) in enumerate(features) if f != '?']
    print ' '.join(features)
