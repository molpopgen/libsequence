#ifndef SEQUENCE_VARIANT_MATRIX_MSFORMAT_HPP__
#define SEQUENCE_VARIANT_MATRIX_MSFORMAT_HPP__

#include <istream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    /// \example ms_to_VariantMatrix.cc

    /*! \brief Create VariantMatrix from "ms"-like input format
     * \param input_stream A model of std::istream
     * \return A VariantMatrix
     * \ingroup variantmatrix
     *
     * See ms_to_VariantMatrix.cc for example.
     */
    template <typename streamtype>
    inline VariantMatrix
    from_msformat(streamtype& input_stream)
    {
        char ch;
        while (!input_stream.eof())
            {
                input_stream >> ch;
                if (ch == ':')
                    break;
            }
        std::string temp;
        std::size_t S;
        input_stream >> S >> temp;
        std::vector<double> pos(S);
        for (std::size_t i = 0; i < S; ++i)
            {
                input_stream >> pos[i];
            }
        std::vector<std::int8_t> data;
        char next_token;
        while (!input_stream.eof()
               && static_cast<char>(input_stream.peek()) != '/')
            {
                input_stream >> next_token >> std::ws;
                data.push_back((next_token == '0') ? 0 : 1);
            }
        input_stream >> std::ws;
        //We now have to transpose the input matrix
        decltype(data) data_t(data.size(), -1);
        std::size_t nsam = data.size() / static_cast<std::size_t>(S);
        std::size_t k = 0;
        for (std::size_t i = 0; i < static_cast<std::size_t>(S); ++i)
            {
                for (std::size_t j = 0; j < nsam; ++j)
                    {
                        data_t[k++] = data[i + j * S];
                    }
            }
        return VariantMatrix(std::move(data_t), std::move(pos));
    }

    template <typename output_stream>
    inline void
    to_msformat(const VariantMatrix& m, output_stream& o)
    /*! \brief Write VariantMatrix in "ms" format.
     * \param m A VariantMatrix
     * \param o A model of std::ostream
     * \ingroup variantmatrix
     *
     * See ms_to_VariantMatrix.cc for example.
     */
    {
        o << "//\nsegsites: " << m.nsites << "\npositions: ";
        for (auto& p : m.positions)
            {
                o << p << ' ';
            }
        o << '\n';
        for (std::size_t i = 0; i < m.nsam; ++i)
            {
                auto col = get_ConstColView(m, i);
                for (auto state : col)
                    {
                        o << static_cast<int>(state);
                    }
                if (i < m.nsam - 1)
                    {
                        o << '\n';
                    }
            }
    }
} // namespace Sequence

#endif
