import glob
import os

MyDir=os.path.join(os.path.dirname(__file__),'example_files')

__all__  = ['MyDir', 'reading_test']

def reading_test(Pattern, Reader, DEBUG=False):
    nError=0
    if DEBUG:
        print('')
    failedTest=[]
    for f in glob.glob(os.path.join(MyDir,Pattern)):
        if os.path.splitext(f)[-1] in ['.py','.pyc'] or f.find('_TMP')>0:
            continue
        try:
            obj = Reader(f)
            s=type(obj).__name__.replace('File','')[:20]
            if DEBUG:
                print('[ OK ] {:30s}\t{:20s} {:20s}'.format(os.path.basename(f)[:30],s))
        except:
            nError += 1
            failedTest=[os.path.basename(f)]
            if DEBUG:
                print('[FAIL] {:30s}\tException occurred'.format(os.path.basename(f)[:30]))
            raise 
    if nError>0:
        print('Reading failed for files: ',failedTest)
        raise Exception('Some tests failed')
