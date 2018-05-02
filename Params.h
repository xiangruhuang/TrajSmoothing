#ifndef PARAMS_H
#define PARAMS_H
class Params{
    public:
        Params(int argc, char** argv){
            parse_cmd_line(argc, argv);
        }

        void exit_with_help(){
            cerr << "./motion_smooth (Options) [filelist path] (output_folder, default to 'motion/smooth')" << endl;
            cerr << "Options:" << endl;
            cerr << "-t: set threshold (200.0)" << endl;
            cerr << "-method: set correspondence method, candidates are {shortest_path}" << endl;
            cerr << "-lr: set learning rate (1e-3)" << endl;
            cerr << "-o: set number of sample points to characterize curves (0)" << endl;
            cerr << "-k: set K for number of nearest neighbors (10)" << endl;
            cerr << "-c_align: set C_align to control alignment (1.0)" << endl;
            cerr << "-c_smooth: set C_smooth to control smoothing (1.0)" << endl;
            cerr << "-c_reg: set C_reg for l2 regularization (1.0)" << endl;
            cerr << "-solver: set solver type: {GD, AM}" << endl;
            cerr << "-tol: set tolerence for termination (0.1)" << endl;
            cerr << "-interval: set length for sampling (0.001)" << endl;
            cerr << "-sigma: set sigma for computing weights (0.05)" << endl;
            cerr << "-radius: set radius for computing correspondence (0.05)" << endl;
            cerr << "-lambda: set lambda for curve smoothing (1.0)" << endl;
            cerr << "-ratio: set ratio for topological reconstruction (0.5)" << endl;
            cerr << "-max_iter: set max number of iteration (30000)" << endl;
            exit(1);
        }

        void parse_cmd_line(int argc, char** argv){
            int i;
            vector<string> args;
            for (int i = 0; i < argc; i++){
                string arg_i(argv[i]);
                args.push_back(arg_i);
            }
            
            for(i=1;i<argc;i++){
                if( args[i][0] != '-' )
                    break;
                if( ++i >= argc )
                    exit_with_help();
                if (args[i-1] == "-t"){
                    threshold = stod(args[i]);
                    continue;
                }
                if (args[i-1] == "-lr"){
                    learning_rate = stod(args[i]);
                    continue;
                }
                if (args[i-1] == "-method"){
                    method = args[i];
                    continue;
                }
                if (args[i-1] == "-tol"){
                    tol = stod(args[i]);
                    continue;
                }
                if (args[i-1] == "-radius"){
                    radius= stod(args[i]);
                    continue;
                }
                if (args[i-1] == "-interval"){
                    interval = stod(args[i]);
                    continue;
                }
                if (args[i-1] == "-o"){
                    ord = stoi(args[i]);
                    continue;
                }
                if (args[i-1] == "-D"){
                    D = stoi(args[i]);
                    continue;
                }
                if (args[i-1] == "-max_iter"){
                    max_iter = stoi(args[i]);
                    continue;
                }

                if (args[i-1] == "-k"){
                    K = stoi(args[i]);
                    continue;
                }
                if (args[i-1] == "-sigma"){
                    sigma = stod(args[i]);
                    continue;
                }
                if (args[i-1] == "-c_align"){
                    C_align = stod(args[i]);
                    continue;
                }
                if (args[i-1] == "-c_smooth"){
                    C_smooth = stod(args[i]);
                    continue;
                }
                if (args[i-1] == "-c_reg"){
                    C_reg = stod(args[i]);
                    continue;
                }
                if (args[i-1] == "-lambda"){
                    lambda = stod(args[i]);
                    continue;
                }

                if (args[i-1] == "-solver"){
                    solver = args[i];
                    continue;
                }

                if (args[i-1] == "-ratio"){
                    ratio = stod(args[i]);
                    continue;
                }
                
                cerr << "unknown option: " << args[i-1] << endl;
                exit(0);
            }

            if(i>=argc)
                exit_with_help();

            filelist_path = args[i];
            if (i + 1 < argc){
                output_folder = args[i+1];
            }
        }

        void dump(){
            cerr << "Parameters: {" << endl;
            cerr << "\tthreshold=" << threshold << endl;
            cerr << "\tfilelist_path=" << filelist_path << endl;
            cerr << "\tord=" << ord << endl;
            cerr << "\tlearning rate=" << learning_rate << endl;
            cerr << "\tK=" << K << endl;
            cerr << "\tD=" << D << endl;
            cerr << "\toutput_folder=" << output_folder << endl;
            cerr << "\tcorrespondence method=" << method << endl;
            cerr << "\tC_reg=" << C_reg << ", C_align=" << C_align << 
                ", C_smooth=" << C_smooth << endl;
            cerr << "\tinterval=" << interval << endl;
            cerr << "\tsolver=" << solver << endl;
            cerr << "\ttol=" << tol << endl;
            cerr << "\tsigma=" << sigma << endl;
            cerr << "\tradius=" << radius << endl;
            cerr << "\tratio=" << ratio << endl;
            cerr << "\tmax_iter=" << max_iter << endl;
            cerr << "}" << endl;
        }

        Float threshold = 100.0;
        int K = 10, ord = 0, D = 2;
        int max_iter = 30000;
        Float learning_rate = 1e-3;
        string filelist_path, output_folder;
        Float C_reg = 1.1, C_align = 1.2, C_smooth = 1.3;
        Float tol = 1e-1;
        Float interval = 0.001;
        Float sigma = 0.005;
        Float radius = 0.05;
        string method;
        string solver;
        bool dump_matchings = true;
        Float lambda = 1.0;
        Float ratio = 0.5;
};
#endif
