#ifndef SIMOL_COMMANDLINE_HPP
#define SIMOL_COMMANDLINE_HPP

#include <getopt.h>

namespace simol
{
    static struct option long_options[] =
    {
      {"input-file", required_argument, 0, 'i'},
      {0, 0, 0, 0}
    };

  class CommandLine
  {
    public:
      CommandLine(int argc, char* argv[])
      {
        int c;

        while (1)
        {
          /* getopt_long stores the option index here. */
          int option_index = 0;

          c = getopt_long (argc, argv, "abc:d:f:i:",
                           long_options, &option_index);

          /* Detect the end of the options. */
          if (c == -1)
            break;

          switch (c)
          {
            case 0:
            /* If this option sets a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
              break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
              printf (" with arg %s", optarg);
            printf ("\n");
            break;

            case 'i':
              inputFileName_ = optarg;
              printf ("option -i with value `%s'\n", optarg);

            case '?':
            /* getopt_long already printed an error message. */
            break;

            default:
              abort ();
          }
        }
      }

    std::string const & inputFileName() const
    { return inputFileName_; }

    private:
      std::string inputFileName_;
  };
}

#endif
