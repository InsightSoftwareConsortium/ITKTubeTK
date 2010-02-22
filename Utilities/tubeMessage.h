
namespace mavPrintLevel
{
  enum { Information, Warning, Error, Debug };
}

template <class T>
void mavPrint( const T& str, int level = 0 )
{
  #ifndef NDEBUG
  if( level == mavPrintLevel::Debug )
    {
    std::cout << "<debug>" << str << "</debug>" << std::endl;
    return;
    }
  #endif
  switch( level )
    {
    case mavPrintLevel::Information:
      std::cout << "<info>" << str << "</info>" << std::endl;
      break;
    case mavPrintLevel::Warning:
      std::cout << "<warning>" << str << "</warning>" << std::endl;
      break;
    case mavPrintLevel::Error:
      std::cout << "<error>" << str << "</error>" << std::endl;
      break;
    default:
      break;
    }
}

template <class T>
void mavPrintInfo( const T& str )
{
  mavPrint( str, mavPrintLevel::Information );
}

template <class T>
void mavPrintWarning( const T& str )
{
  mavPrint( str, mavPrintLevel::Warning );
}

template <class T>
void mavPrintError( const T& str )
{
  mavPrint( str, mavPrintLevel::Error );
}

template <class T>
void mavPrintDebug( const T& str )
{
  mavPrint( str, mavPrintLevel::Debug );
}

template <class T>
void mavAssert( bool assertion, const T& str, int level=mavPrintLevel::Error )
{
  if( level == mavPrintLevel::Debug && !assertion )
    {
    #ifndef NDEBUG
      mavPrintDebug( str );
      mav::Exception e;
      throw e;
    #endif
    return;
    }
  if( level == mavPrintLevel::Information && !assertion )
    {
    mavPrintInfo( str );
    #ifndef NDEBUG
      mav::Exception e;
      throw e;
    #endif
    return;
    }
  if( !assertion ) // Warning or Error
    {
    mavPrint( std::string( "mavAssertion Failed: " ) + mavToString( str ),
              level );
    mav::Exception e;
    throw e;
    }
}

template <class T>
void mavAssertInfo( bool assertion, const T& str )
{
  mavAssert( assertion, str, mavPrintLevel::Information );
}

template <class T>
void mavAssertWarning( bool assertion, const T& str )
{
  mavAssert( assertion, str, mavPrintLevel::Warning );
}

template <class T>
void mavAssertError( bool assertion, const T& str )
{
  mavAssert( assertion, str, mavPrintLevel::Error );
}

template <class T>
void mavAssertDebug( bool assertion, const T& str )
{
  mavAssert( assertion, str, mavPrintLevel::Debug );
}

