/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#ifndef __BITPIT_BINARY_STREAM_HPP__
#define __BITPIT_BINARY_STREAM_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdexcept>

// Bitpit
// none

// Others
// none

// ========================================================================== //
// TYPES DEFINITIONS                                                          //
// ========================================================================== //
// none

// Forward declarations ------------------------------------------------- //
namespace bitpit{
class IBinaryStream;
class OBinaryStream;
};


// Function prototypes -------------------------------------------------- //
template<typename T>
bitpit::IBinaryStream& operator>>(                                                  // Stream operator for class IBinaryStream
    bitpit::IBinaryStream                     &istm,                                // (input) input stream
    T                               &val                                  // (input) value to be streamed
);

template<typename T>
bitpit::OBinaryStream& operator<<(                                                  // Stream operator for class OBinaryStream
    bitpit::OBinaryStream                     &ostm,                                // (input) output stream
    const T                         &val                                  // (input) value to be streamed
);
template<>
bitpit::OBinaryStream& operator<<(                                                  // Explicit specialization of input stream operator for std::string
    bitpit::OBinaryStream                     &ostm,                                // (input) output stream
    const std::string               &val                                  // (input) string to be streamed
);
bitpit::OBinaryStream& operator<<(                                                  // Stream operator for class OBinaryStream
    bitpit::OBinaryStream                     &ostm,                                // (input) output stream
    const char                      *val                                  // (input) char array to be streamed
);

namespace bitpit{

// Class IBinaryStream ---------------------------------------------------- //
class IBinaryStream {

    // Member(s) ======================================================== //
    private:

    std::vector<char>               buffer;                               // stream buffer
    size_t                          current_pos;                          // Cursor position

    // Constructor(s) =================================================== //
    public:

    IBinaryStream(                                                          // Default constructor (empty stream)
        void                                                              // (input) none
    );
    IBinaryStream(                                                          // Custom constructor #1 (empty stream with known capacity);
        size_t                       capacity                               // (input) buffer capacity
    );
    IBinaryStream(                                                          // Custom constructor #2 (stream pointing to memory location)
        const char*                  buf_,                                // (input) pointer to memory location
        size_t                       capacity                             // (input) buffer capacity
    );
    IBinaryStream(                                                          // Custom constructor #3 (stream initialized from std::vector<char>)
        const std::vector<char>          &vec                                  // (input) vector used for initialization
    );

    // Destructor(s) ==================================================== //
    // default

    // Public method(s) ================================================= //
    public:
    void setCapacity(                                                     // Set the capacity of the stream
        size_t                       capacity                             // (input) new capacity (in bytes) of stream
    );
    size_t capacity(                                                  // Capacity of the stream
        void
    ) const;
    void open(                                                            // Open input stream from memory location
        const char                  *mem,                                 // (input) pointer to memory location
        size_t                       capacity                             // (input) capacity (in bytes) of memory chunk
    );
    bool eof(                                                             // Flag for eof
        void                                                              // (input) none
    ) const;
    std::ifstream::pos_type tellg(                                             // Returns current position of cursor in the buffer
        void                                                              // (input) none
    ) const;
    bool seekg (                                                          // Set cursor position in the current buffer
        size_t                       pos                                  // (input) position
    );
    bool seekg (                                                          // Set cursor position in the current buffer
        std::streamoff               offset,                              // (input) offset with respect to the specified direction
        std::ios_base::seekdir       way                                  // (input) offset direction
    );
    const std::vector<char>& data(                                 // Returns reference to buffer
        void                                                              // (input) none
    ) { return(buffer); }
    char* rawData(                                                     // Returns pointer to buffer
        void                                                              // (input) none
    ) { return( buffer.data() ); }

    // Private methods(s) =============================================== //
    private:

    template<typename T>
    void read(                                                            // Read data from memory location pointed by t and store into stream buffer
        T                           &t                                    // (input) data to be imported in the stream buffer
    );
    void read(                                                            // Explicit template specialization of IBinaryStream::read for vector<char>
        std::vector<char>           &vec                                  // (input) source vector
    );
    void read(                                                                // Read data from memory location pointed by p and store into stream buffer
        char                        *p,                                   // (input) pointer to memory location
        size_t                       size                                 // (input) size (in bytes) of data to be read
    );

    // Friendships ====================================================== //
    template< typename T >
    friend IBinaryStream& (::operator >>) (IBinaryStream&, T& );
};

// Class OBinaryStream ---------------------------------------------------- //
class OBinaryStream {

    // Member(s) ======================================================== //
    private:

    size_t                           current_pos;                         // Cursor current position
    std::vector<char>                buffer;                              // Buffer

    // Constructor(s) =================================================== //
    public:

    OBinaryStream(                                                          // Default constructor (create empty object)
        void                                                              // (input) none
    );
    OBinaryStream(                                                          // Custom constructor #1 (create an empty object with buffer of specified capacity)
        size_t                       capacity                               // (input) none
    );

    // Destructor(s) ==================================================== //
    // none

    // Public method(s) ================================================= //
    void setCapacity(                                                     // Set the capacity of the stream
        size_t                       capacity                             // (input) new capacity (in bytes) of stream
    );
    size_t capacity(                                                      // Capacity of the stream
        void
    ) const;
	void open(                                                            // Open output stream
        size_t                       capacity                             // (input) stream capacity
    );
    bool eof(                                                             // Flag for eof
        void                                                              // (input) none
    ) const;
    std::ifstream::pos_type tellg(                                             // Returns current position of cursor in the buffer
        void                                                              // (input) none
    ) const;
    bool seekg (                                                          // Set cursor position in the current buffer
        size_t                       pos                                  // (input) position
    );
    bool seekg (                                                          // Set cursor position in the current buffer
        std::streamoff               offset,                              // (input) offset with respect to the specified direction
        std::ios_base::seekdir       way                                  // (input) offset direction
    );
    void squeeze(                                                         // Squeeze the stream to fit the data
        void                                                              // (input) none
    );
    const std::vector<char>& data(                                        // Returns reference to buffer
        void                                                              // (input) none
    ) { return(buffer); }
    char* rawData(                                                     // Returns pointer to buffer
        void                                                              // (input) none
    ) { return( buffer.data() ); }

    // Private method(s) ================================================ //
    private:

    template<typename T>
    void write(                                                           // Write data to internal buffer
        const T                     &t                                    // (input) data to be written in the internal buffer
    );
    void write(                                                           // Write char array to internal buffer
        const char                  *p,                                   // (input) pointer to char array
        size_t                       size                                 // (input) size of data chunk to be written in the internal buffer
    );
    void write(                                                           // Write vector of char to internal buffer
        const std::vector<char>     &vec                                  // (input) vector of char to be written in the internal buffer
    );

    // Friendship(s) ==================================================== //
    template<typename T>
    friend OBinaryStream& (::operator<<) ( OBinaryStream&, const T& );
};

}

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
# include "binary_stream.tpp"

#endif
