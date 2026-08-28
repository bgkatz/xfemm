#ifndef FILETOKENIZER_H
#define FILETOKENIZER_H

#include <cstdio>
#include <cstdlib>
#include <vector>

namespace femm {

/**
 * @brief Reads a whole text file into memory and hands out
 * whitespace-separated numeric tokens.
 *
 * This is dramatically faster than per-token fscanf (which locks the
 * stream and re-parses the format string on every call) for the large
 * mesh files (.node/.ele/.edge) written by the mesher.
 */
class FileTokenizer
{
public:
    explicit FileTokenizer(const char *path)
        : pos(nullptr)
        , fileEnd(nullptr)
    {
        FILE *fp = fopen(path, "rb");
        if (fp == nullptr)
            return;
        fseek(fp, 0, SEEK_END);
        long size = ftell(fp);
        fseek(fp, 0, SEEK_SET);
        if (size >= 0)
        {
            buf.resize((size_t)size + 1);
            size_t nread = fread(buf.data(), 1, (size_t)size, fp);
            buf[nread] = '\0';
            pos = buf.data();
            fileEnd = buf.data() + nread;
        }
        fclose(fp);
    }

    /// \brief returns true if the file was read into memory successfully
    bool isOpen() const
    {
        return pos != nullptr;
    }

    /// \brief read the next whitespace-delimited integer; returns false at end of file
    bool nextInt(int &v)
    {
        skipSpace();
        if (pos >= fileEnd)
            return false;
        char *end;
        v = (int)strtol(pos, &end, 10);
        if (end == pos)
            return false;
        pos = end;
        return true;
    }

    /// \brief read the next whitespace-delimited floating point value; returns false at end of file
    bool nextDouble(double &v)
    {
        skipSpace();
        if (pos >= fileEnd)
            return false;
        char *end;
        v = strtod(pos, &end);
        if (end == pos)
            return false;
        pos = end;
        return true;
    }

    /// \brief skip the remainder of the current line (e.g. unused header fields)
    void skipLine()
    {
        while (pos < fileEnd && *pos != '\n')
            ++pos;
        if (pos < fileEnd)
            ++pos;
    }

private:
    void skipSpace()
    {
        while (pos < fileEnd &&
               (*pos == ' ' || *pos == '\t' || *pos == '\r' || *pos == '\n'))
            ++pos;
    }

    std::vector<char> buf;
    char *pos;
    const char *fileEnd;
};

} // namespace femm

#endif
