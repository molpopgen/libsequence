// Code for the -*- C++ -*- namespace Sequence::ClustalW<T>

/*

Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

Remove the brackets to email me.

This file is part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.

*/

/*! \file Clustalw.tcc
  @brief code for Clustalw.hpp
*/
#include <map>
#include <Sequence/AlignStream.hpp>
#include <iterator>
#include <algorithm>
#include <cassert>

namespace Sequence
{
    template <typename T>
    std::istream &
    ClustalW<T>::read(std::istream &s)
    /*!
    Calls to Sequence::operator>> into objects of type ClustalW<T>
    results in a call to this function, which reads the alignment in
    from the stream.
  */
    {
        std::string clustalw;
        char ch;
        std::map<std::string, std::string> seqs;
        std::map<std::string, int> order;
        std::size_t nseqs = 0;
        s >> clustalw >> std::ws;
        if (clustalw != "CLUSTAL")
            {
                throw std::runtime_error("Sequence::ClustalW::read() : input "
                                         "stream does not appear to be in "
                                         "CLUSTALW format");
            }
        else
            {
                ReadThroughLine(s);
            }
        std::string temp, temp2;
        while (!s.eof())
            {
                s.get(ch);
                bool putback = 0;
                if (ch == '\n' || ch == ' ' || ch == '*')
                    {
                        s.putback(ch);
                        ReadThroughLine(s);
                        putback = 1;
                    }
                else
                    {
                        if (!putback)
                            s.putback(ch);
                        s >> temp;
                        auto iter = seqs.find(temp);
                        if (iter != seqs.end())
                            {
                                s >> temp2 >> std::ws;
                                seqs[(*iter).first] += temp2;
                            }
                        else
                            {
                                s >> temp2 >> std::ws;
                                seqs[temp] = temp2;
                                order[temp] = nseqs++;
                            }
                    }
                s >> std::ws;
            }

        typename std::vector<T> _data;
        for (int i = 0; i < nseqs; ++i)
            {
                auto iter = seqs.begin(), iter_end = seqs.end();
                bool found = 0;
                while (iter != iter_end)
                    {
                        if (order[(*iter).first] == i)
                            {
                                T t;
                                t.name = iter->first;
                                t.seq = std::move(iter->second);
                                _data.emplace_back(std::move(t));
                                //_data[i].name = std::move(iter->first);
                                //_data[i].seq = std::move(iter->second);
                                iter = iter_end;
                                found = 1;
                            }
                        if (!found)
                            ++iter;
                    }
            }
        if (_data.size() != nseqs)
            {
                throw std::runtime_error("fatal error converting input data");
            }
        this->assign(std::move(_data));
        return s;
    }

    template <typename T>
    std::ostream &
    ClustalW<T>::print(std::ostream &s) const
    {
        typename ClustalW<T>::const_iterator i = this->begin(),
                                             j = this->end();
        unsigned k = 0, len = unsigned(i->seq.length());
        s << "CLUSTAL W"
          << "\n\n";
        while (k < len)
            {
                unsigned offset = (k + 60 < len) ? k + 60 : k + (len - k);
                for (i = this->begin(); i < j; ++i)
                    {
                        s << i->name << '\t';
                        std::copy(i->seq.begin()
                                      + std::string::difference_type(k),
                                  i->seq.begin()
                                      + std::string::difference_type(offset),
                                  std::ostream_iterator<char>(s, ""));
                        s << '\n';
                    }
                s << '\n';
                k = offset;
            }
        return s;
    }

    template <typename T>
    std::istream &
    ClustalW<T>::ReadThroughLine(std::istream &s)
    {
        char ch;
        while (s.get(ch))
            {
                if (ch == '\n')
                    return s;
            }
        return s;
    }
}
