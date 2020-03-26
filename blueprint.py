#!/usr/bin/env python3
"""
Generate a blueprint file from a pdb file.
Optionally including remodelling sites.
"""

amino_acids = {
        "ALA":"A","GLY":"G","ILE":"I",
        "LEU":"L","PRO":"P","VAL":"V",
        "PHE":"F","TRP":"W","TYR":"Y",
        "ASP":"D","GLU":"E","ARG":"R",
        "HIS":"H","LYS":"K","SER":"S",
        "THR":"T","CYS":"C","MET":"M",
        "ASN":"N","GLN":"Q"}

def read_pdb_sequence(pdb_file):
    """
    Read the one letter sequence from a pdb file
    """
# PDB format:
#ATOM      1  N   ALA P   1     -65.001  34.632  12.643  1.00  0.00      PA1A
    with open(pdb_file) as fh:
        resnum = 0
        for line in fh:
            if not line.startswith("ATOM"): continue
            resname = amino_acids[line[17:20]] # pdb pos
            new_resnum = int(line[22:26]) # pdb pos
            if resnum == new_resnum: continue
            resnum = new_resnum
            yield resnum,resname

def make_chunks(*iterables, nchunks=1):
    output = []
    for iterable in iterables:
        iterable = list(iterable)
        chunkout = []
        length, rest = divmod(len(iterable), nchunks)
        for i in range(nchunks - 1):
            chunkout.append(iterable[0+i*length:length+i*length])
        chunkout.append(iterable[-length-rest+1:])
        output.append(chunkout)
    output = list(zip(*output))
    return output

class Chain():
    """
    Collection of the data that describes the chain
    """
    def __init__(self, resnums, resnames, ss=None, bp=None):
        self.bp = bp
        self.resnums = resnums
        self.resnames = resnames
        if ss is None:
            nres = len(self.resnums)
            assert nres > 0
            self.ss = ["."] * nres
        else:
            self.ss = ss

    def splice(self, tosplice, after):
        """
        splice in a piece of chain
        """
        after = self.bp.index[after]
        self.resnums = insert_list(self.resnums, after, tosplice.resnums)
        self.resnames = insert_list(self.resnames, after, tosplice.resnames)
        self.ss = insert_list(self.ss, after, tosplice.ss)


def insert_list(acceptor, after, donor):
    left = acceptor[:after]
    right = acceptor[after:]
    return left + donor + right

class Blueprint():
    """

    """
    def __init__(self, resnums, resnames, ss=None, nchains=1):
        if nchains > 1:
            self.chains = []
            for i, [c_resnums, c_resnames] in enumerate(make_chunks(resnums, resnames, nchunks=nchains)):
                # I see ss is missing here...
                self.chains.append(Chain(c_resnums, c_resnames, bp=self))
        else:
            self.chains = [Chain(resnums, resnames, ss, bp=self)]
        self.make_index()

    def make_index(self):
        # This will fail when there is N-terminal extension
        index = dict()
        for chain in self.chains:
            index.update({ri: idx for ri, idx in zip(chain.resnums, list(range(len(chain.resnums)))) if ri != 0})
        self.index = index

    @classmethod
    def from_pdb(cls, filename, *args, **kwargs):
        num_seq = read_pdb_sequence(filename)
        resnums = []
        resnames = []
        for ri, rn in num_seq:
            resnums.append(ri)
            resnames.append(rn)
        return cls(resnums, resnames, *args, **kwargs)

    @classmethod
    def from_file(cls):
        raise NotImplementedError

    def write(self, outfile=f"remodel.bp"):
        """
        Write the blueprint to file
        """
        sequential = 1 # Rosetta numbering is sequential starting from 1
        with open(outfile, 'w') as fh:
            for chain in self.chains:
                for num, name, ss in zip(chain.resnums, chain.resnames, chain.ss):
                    # This is getting super confusing...
                    if num != 0:
                        # This will totally break for N-terminal extension
                        num = sequential
                        sequential += 1
                    fh.write(" ".join([str(num), name, ss])+"\n")

        print(f"Wrote {outfile}")

    def parse_extend(self, extendstr):
        """
        Parse an extend command
        Form of:
            after:length:ss
        """
        for cmd in extendstr.split(","):
            cmd = cmd.split(":")
            assert len(cmd) == 3
            after, length, ss = cmd
            after = int(after)
            length = int(length)
            ss = str(ss)
            self.extend(after, length, ss)

    def extend(self, after, length, ss):
        """
        Low level extension
        Make an extension
        """
        # Don't forget, N-terminal should only be like this
        # In the first chain
        if after == "N":
            raise NotImplementedError
        # C-terminal only in the last
        if after == "C":
            raise NotImplementedError
        # Lets start with in the middle because i need that
        # Lets start with unknown sequence but known ss
        # format is:
        # ri rn ss
        # 0 x H
        # number always 0
        # name always x (or i suppose this can be specified..)
        # ss has to be specified
        resnums = [0] * length
        resnames = ["x"] * length
        l_ss = [ss] * length
        tosplice = Chain(resnums, resnames, l_ss)
        for chain in self.chains:
            chain.splice(tosplice, after=after)
        self.make_index()

        # Will do termini later

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-pdb", default=None)
    parser.add_argument("-bp", default=None)
    parser.add_argument("-extend", type=str, default=None)
    parser.add_argument("-nchains", type=int, default=1)
    parser.add_argument("-outfile", default="remodel.bp")

    args = parser.parse_args()

    # Initialize a blueprint
    blueprint = None
    if args.pdb is not None:
        blueprint = Blueprint.from_pdb(
            args.pdb,
            nchains=args.nchains)
    elif args.bp is not None:
        blueprint = Blueprint.from_file(
            args.bp,
            nchains=args.nchains)
    elif blueprint is None:
        print("No blueprint initialized")
        parser.print_help()
        exit()

    # Do optional modifications
    if args.extend is not None:
        blueprint.parse_extend(args.extend)

    blueprint.write(args.outfile)
