def read_primary_header(fn):
    import fitsio
    '''
    Reads the FITS primary header (HDU 0) from the given filename.
    This is just a faster version of fitsio.read_header(fn).
    '''
    if fn.endswith('.gz'):
        return fitsio.read_header(fn)

    # Weirdly, this can be MUCH faster than letting fitsio do it...
    hdr = fitsio.FITSHDR()
    foundEnd = False
    ff = open(fn, 'rb')
    h = b''
    while True:
        hnew = ff.read(32768)
        if len(hnew) == 0:
            # EOF
            ff.close()
            raise RuntimeError('Reached end-of-file in "%s" before finding end of FITS header.' % fn)
        h = h + hnew
        while True:
            line = h[:80]
            h = h[80:]
            #print('Header line "%s"' % line)
            # HACK -- fitsio apparently can't handle CONTINUE.
            # It also has issues with slightly malformed cards, like
            # KEYWORD  =      / no value
            if line[:8] != b'CONTINUE':
                try:
                    hdr.add_record(line.decode())
                except OSError as err:
                    print('Warning: failed to parse FITS header line: ' +
                          ('"%s"; error "%s"; skipped' % (line.strip(), str(err))))
                          
            if line == (b'END' + b' '*77):
                foundEnd = True
                break
            if len(h) < 80:
                break
        if foundEnd:
            break
    ff.close()
    return hdr
