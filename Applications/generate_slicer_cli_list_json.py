import os, sys, json

if __name__ == "__main__":

  exclude_cli_list = [
    'ImageMath',
    'TubeMath',
    'TreeMath'
    ]

  appdir = os.path.dirname(os.path.realpath(__file__))

  cli_list = [f
              for f in os.listdir(appdir)
              if (os.path.isdir(f) and
                  f not in exclude_cli_list and
                  os.path.isfile(os.path.join(appdir, f, f + '.xml'))
                  )]

  slicer_cli_list_json = {}

  for cli in cli_list:

    if os.path.isfile(os.path.join(appdir, cli, cli + '.py')):

      slicer_cli_list_json[cli] = {'type': 'python'}

    else:

      slicer_cli_list_json[cli] = {'type': 'cxx'}

  with open(os.path.join(appdir, 'slicer_cli_list.json'), 'w') as f:
    json.dump(slicer_cli_list_json, f, indent=1, sort_keys=True)

