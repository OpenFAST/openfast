from pathlib import Path

from bs4 import BeautifulSoup

with open('configsurations.xml') as fp:
    configs = soup = BeautifulSoup(fp, 'xml')

for path in Path('.').rglob('*.vfproj'):
    print(path)
    with open(path) as fp:
        soup = BeautifulSoup(fp, 'xml')
        for cfg in soup.findAll('Configuration'):
            configs[cfg['Name']] = cfg
    break

print(configs)