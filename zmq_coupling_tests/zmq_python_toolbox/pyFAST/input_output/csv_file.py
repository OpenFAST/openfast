import os

from .file import File, WrongFormatError
import pandas as pd

class CSVFile(File):
    """ 
    Read/write a CSV file. 

    Main methods
    ------------
      read, write, toDataFrame

    Examples
    --------

        # Read a csv file and convert it to a pandas dataframe
        f = CSVFile('test.csv')
        df = f.toDataFrame()

    """

    @staticmethod
    def defaultExtensions():
        return ['.csv','.txt']

    @staticmethod
    def formatName():
        return 'CSV file'

    def __init__(self, filename=None, sep=None, colNames=None, commentChar=None, commentLines=None,\
                       colNamesLine=None, detectColumnNames=True, header=None, **kwargs):
        colNames     = [] if colNames is None else colNames
        commentLines = [] if commentLines is None else commentLines
        self.sep          = sep
        self.colNames     = colNames
        self.commentChar  = commentChar
        self.commentLines = commentLines
        self.colNamesLine = colNamesLine
        self.detectColumnNames = detectColumnNames
        self.data=[]
        if header is None:
            self.header=[]
        else:
            if not hasattr(header, '__len__'):
                self.header=[header]
            else:
                self.header=header
        self.nHeader=0
        if (len(self.commentLines)>0) and (self.commentChar is not None):
            raise Exception('Provide either `commentChar` or `commentLines` for CSV file types')
        if (len(self.colNames)>0) and (self.colNamesLine is not None):
            raise Exception('Provide either `colNames` or `colNamesLine` for CSV file types')
        super(CSVFile, self).__init__(filename=filename,**kwargs)

    def _read(self):
        COMMENT_CHAR=['#','!',';']
        # --- Detecting encoding
        # NOTE: done by parent class method
        
        # --- Subfunctions
        def readFirstLines(nLines):
            lines=[]
            with open(self.filename, 'r', encoding=self.encoding, errors="surrogateescape") as fid:
                for i, line in enumerate(fid):
                    lines.append(line.strip())
                    if i==nLines:
                        break
            return lines

        def readline(iLine):
            with open(self.filename,'r',encoding=self.encoding) as f:
                for i, line in enumerate(f):
                    if i==iLine:
                        return line.strip()
                    elif i>iLine:
                        break
        def split(s):
            if s is None:
                return []
            if self.sep==r'\s+':
                return s.strip().split()
            else:
                return [c.strip() for c in s.strip().split(self.sep)]
        def strIsFloat(s):
            try:
                float(s)
                return True
            except:
                return False
        # --- Safety
        if self.sep=='' or self.sep==' ':
            self.sep=r'\s+'

        iStartLine=0
        
        # --- Exclude some files from the CSV reader ---
        line=readline(iStartLine)
        words=line.split()
        if len(words)>1:
            try:
                int(words[0])
                word0int = True
            except:
                word0int = False
            if word0int and words[1].isalpha():
                raise WrongFormatError('Input File {}: '.format(self.filename) + 'is not likely a CSV file' )
                
        # --- Headers (i.e. comments)
        # TODO: read few headers lines instead of multiple read below..

        self.header = []
        if len(self.commentLines)>0:
            # We read the lines
            with open(self.filename,'r',encoding=self.encoding) as f:
                for i in range(max(self.commentLines)+1):
                    l = f.readline()
                    if i in self.commentLines:
                        self.header.append(l.strip())
        elif self.commentChar is not None:
            # we detect the comments lines that start with comment char
            with open(self.filename,'r',encoding=self.encoding) as f:
                n=0
                while n<100:
                    l = f.readline().strip()
                    if (not l) or (l+'_dummy')[0] != self.commentChar[0]:
                        break
                    self.header.append(l.strip())
                    n+=1
            self.commentLines=list(range(len(self.header)))
        else:
            # We still believe that some characters are comments
            line=readline(iStartLine)
            line=str(line).strip()
            if len(line)>0 and line[0] in COMMENT_CHAR:
                self.commentChar=line[0]
                # Nasty copy paste from above
                with open(self.filename,'r',encoding=self.encoding) as f:
                    n=0
                    while n<100:
                        l = f.readline().strip()
                        if (not l) or (l+'_dummy')[0] != self.commentChar[0]:
                            break
                        self.header.append(l.strip())
                        n+=1

        iStartLine = len(self.header)

        # --- File separator 
        if self.sep is None:
            # Detecting separator by reading first lines of the file
            try:
                with open(self.filename,'r',encoding=self.encoding) as f:
                    dummy=[next(f).strip() for x in range(iStartLine)]
                    head=[next(f).strip() for x in range(2)]
                # comma, semi columns or tab
                if head[1].find(',')>0:
                    self.sep=','
                elif head[1].find(';')>0:
                    self.sep=';'
                elif head[1].find('\t')>0:
                    self.sep='\t'
                else:
                    self.sep=r'\s+'
            except:
                # most likely an empty file
                pass

        # --- ColumnNames
        if self.colNamesLine is not None:
            if self.colNamesLine<0:
                # The column names are hidden somwhere in the header 
                line=readline(iStartLine+self.colNamesLine).strip()
                # Removing comment if present (should be present..)
                if self.commentChar is not None:
                    if line.find(self.commentChar)==0:
                        line=line[len(self.commentChar):].strip()
                self.colNames = split(line)
            else:
                line=readline(self.colNamesLine)
                self.colNames=split(line)
                iStartLine = max(iStartLine,self.colNamesLine+1)
        elif len(self.colNames)>0:
            pass
        elif not self.detectColumnNames:
            pass
        else:
            # Looking at first line of data, if mainly floats -> it's not the column names
            colNames = split(readline(iStartLine))
            nFloat = sum([strIsFloat(s) for s in colNames])
            if nFloat ==0 or (len(colNames)>2 and nFloat <= len(colNames)/2):
                # We assume that the line contains the column names
                self.colNames=colNames
                self.colNamesLine = iStartLine
                iStartLine = iStartLine+1
                #  --- Now, maybe the user has put some units below
                first_line = readline(iStartLine)
                #print('>>> first line',first_line)
                first_cols = split(first_line)
                nFloat = sum([strIsFloat(s) for s in first_cols])
                nPa    = first_line.count('(')+first_line.count('[')
                #if nFloat == 0 or nPa>len(self.colNames)/2:
                if nPa>len(self.colNames)/2:
                    # that's definitely some units
                    if len(first_cols)==len(self.colNames):
                        self.colNames=[c.strip()+'_'+u.strip() for c,u in zip(self.colNames, first_cols)]
                    iStartLine = iStartLine+1
            elif len(self.header)>0:
                # Maybe the columns names are in the header
                if self.sep is not None:
                    first_line = readline(iStartLine)
                    first_cols = split(first_line)
                    #print('CommentChar:',self.commentChar)
                    #print('First line:',first_line)
                    #print('First col :',first_cols)
                    for l in self.header:
                        if self.commentChar is not None:
                            if len(self.commentChar)>0:
                                l=l[len(self.commentChar):]
                        cols=split(l)
                        nFloat = sum([strIsFloat(s) for s in cols])
                        if len(cols)==len(first_cols) and nFloat <= len(colNames)-1:
                            self.colNames = cols
                            break
        # --- Reading data
        skiprows = list(range(iStartLine))
        if (self.colNamesLine is not None):
            skiprows.append(self.colNamesLine)
        if (self.commentLines is not None) and len(self.commentLines)>0:
            skiprows = skiprows + self.commentLines
        skiprows =list(sorted(set(skiprows)))
        if self.sep is not None:
            if self.sep=='\t':
                self.sep=r'\s+'
        #print(skiprows)
        try:
#             self.data = pd.read_csv(self.filename,sep=self.sep,skiprows=skiprows,header=None,comment=self.commentChar,encoding=self.encoding)
            with open(self.filename,'r',encoding=self.encoding) as f:
                self.data = pd.read_csv(f,sep=self.sep,skiprows=skiprows,header=None,comment=self.commentChar)
        except pd.errors.ParserError as e:
            raise WrongFormatError('CSV File {}: '.format(self.filename)+e.args[0])

        if (len(self.colNames)==0) or (len(self.colNames)!=len(self.data.columns)):
            self.colNames=['C{}'.format(i) for i in range(len(self.data.columns))]
        self.data.columns = self.colNames;
        self.data.rename(columns=lambda x: x.strip(),inplace=True)

    def _write(self):
        # --- Safety
        if self.sep==r'\s+' or self.sep=='':
            self.sep='\t'
        # Write
        if len(self.header)>0:
            with open(self.filename, 'w', encoding='utf-8') as f:
                f.write('\n'.join(self.header)+'\n')
            with open(self.filename, 'a', encoding='utf-8') as f:
                try:
                    self.data.to_csv(f,   sep=self.sep,     index=False,header=False, line_terminator='\n')
                except TypeError:
                    print('[WARN] CSVFile: Pandas failed, likely encoding error. Attempting a quick and dirty fix.')
                    s=''
                    vals=self.data.values
                    for l in vals:
                        sLine=(self.sep).join([str(v) for v in l])
                        s+=sLine+'\n'
                    f.write(s)
        else:
            self.data.to_csv(self.filename,sep=self.sep,index=False)

    def __repr__(self):
        s = 'CSVFile: {}\n'.format(self.filename)
        s += 'sep=`{}` commentChar=`{}`\ncolNamesLine={}'.format(self.sep,self.commentChar,self.colNamesLine)
        s += ', encoding={}'.format(self.encoding)+'\n'
        s += 'commentLines={}'.format(self.commentLines)+'\n'
        s += 'colNames={}'.format(self.colNames)
        s += '\n'
        if len(self.header)>0:
            s += 'header:\n'+ '\n'.join(self.header)+'\n'
        if len(self.data)>0:
            s += 'size: {}x{}'.format(len(self.data),len(self.data.columns))
        return s

    def _toDataFrame(self):
        return self.data

