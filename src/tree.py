'''Module for handling suffix trees'''
from __future__ import annotations
from collections import deque
from typing import Iterable
from dataclasses import dataclass


@dataclass
class Knæ:
    '''
    Class for storing nodes. \n
    Here Knæ is the node, parent is the Knæ one step toward the root, children is a dictionary of Knæ one step away from the root (where the key in the dictionary represents first letter in the ben of the child) and ben a tuple refering to some interval in a string.
    '''
    parent: Knæ | None
    ben: tuple[int, int] | None
    children: dict[str, Knæ] | int


class SuffixTree:
    '''
    Class for storing and manipulating suffix trees. 
    '''
    root: Knæ
    x: str

    def __init__(self, x: str) -> None:
        '''
        The class SuffixTree is initialized by being created from a string x. Knæ are created so that Knæ.ben is a tuple including first index and excluding last index. This is the naive implementation. \n
        Creating suffix tree T for string x:
        >>> T=SuffixTree(x)
        '''
        self.x = x + '$'
        self.root = Knæ(None, (0, 0), {})
        n = len(self.x)
        current = self.root
        for i, _ in enumerate(x):
            lam = 0
            j = 0
            while i < n:
                if lam == current.ben[1]-current.ben[0] and current.ben[1] != n-1:
                    if self.x[i] in current.children:
                        current = current.children[self.x[i]]
                        lam = 0
                    else:
                        current.children[self.x[i]] = Knæ(
                            current, (i, n-1), i-j)
                        current = self.root
                        break
                if x[current.ben[0]+lam] == self.x[i]:
                    i += 1
                    j += 1
                    lam += 1
                else:
                    current.parent.children[self.x[current.ben[0]]] = Knæ(
                        current.parent, (current.ben[0], current.ben[0]+lam), {x[current.ben[0]+lam]: current})
                    current.parent = current.parent.children[self.x[current.ben[0]]]
                    current.ben = (current.ben[0]+lam, current.ben[1])
                    current.parent.children[self.x[i]] = Knæ(
                        current.parent, (i, n-1), i-j)
                    current = self.root
                    break

    def bft(self, knæ: Knæ | None = None) -> Iterable[int]:
        '''
        Generator for breadth-first-traversal of node (class: Knæ) yielding all integers (indices) downstream of that node. If no node (Knæ) is provided, the function simply does a breadth-first-traversal of the entire tree.
        '''
        kø = deque([])
        if knæ is None:
            kø.append(self.root)
        else:
            kø.append(knæ)
        while kø:
            element = kø.popleft()
            match element:
                case Knæ(_, _, child):
                    if type(child) == int:
                        kø.append(child)
                    else:
                        kø.extend([(child[key]) for key in child])
                case int():
                    yield element
                case _:
                    pass

    def search_for_pattern(self, p: str) -> list[int]:
        '''
        Searches for exact pattern matches of p in self.x of the SuffixTree.\n
        Returns empty list, [], if pattern is the empty string or does not occur in x. E.g.: T=SuffixTree("ABAAB")
        >>> T.search_for_pattern("AB")
        Returns [0,3]
        >>> T.search_for_pattern("")
        Returns []
        >>> T.search_for_pattern("BB")
        Returns []
        '''
        current = self.root
        P = len(p)
        if P == 0:
            return []
        i = 0
        lam = 0
        while i < P:
            if type(current.children) == int:
                pass
            elif lam == current.ben[1]-current.ben[0] and p[i] in current.children:
                current = current.children[p[i]]
                lam = 0
            if self.x[current.ben[0]+lam] == p[i]:
                i += 1
                lam += 1
            else:
                return []
        return list(self.bft(current))
    
    
    def search_approx(self, p, edits):
        good_list=[]
        def search_approx_pattern(node:Knæ, p:str, edits:int,k:int, cigar:list, j:int, i:int)-> list[list,str]:
            print(edits, k, j, i)
            print("".join(cigar))
            u = node.ben[0]
            v = node.ben[1]
            P = len(p)
            if edits<0:
                print("edits<0")
                return
            if j>=P:
                good_list.append([list(self.bft(node)),"".join(cigar)])
                print("HIT!")
                return
            if self.x[u+i]=="$" and type(node.children) == int:
                print("reached end of branch")
                return
            if i == v-u:
                print("This is a node")
                for n in node.children:
                    search_approx_pattern(node.children[n],p, edits, k, cigar, j, 0)
                return

            #Match/substitution:
            if self.x[u+i]==p[j]:
                cigar.append("M")
                search_approx_pattern(node,p, edits,k+1, cigar, j+1, i+1)
            else:
                cigar.append("M")
                search_approx_pattern(node,p, edits-1,k+1, cigar, j+1, i+1)
            #Insertion:
            cigar[k]="I"
            search_approx_pattern(node,p, edits-1,k+1, cigar, j, i+1)

            #Deletion:
            cigar[k]="D"
            search_approx_pattern(node,p, edits-1,k+1, cigar, j, i)
            return
        search_approx_pattern(self.root,p,edits,0,[],0,0)
        return good_list