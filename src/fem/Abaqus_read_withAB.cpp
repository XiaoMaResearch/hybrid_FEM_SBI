//
//  Abaqus_read_withAB.cpp
//  hybrid_fem_bie
//
//  Created by Max on 3/23/18.
//
//

#include "Abaqus_read_withAB.hpp"

void Abaqus_read_withAB(std::string &fname, MatrixXd &Nodes, MatrixXi_rm &Element, std::vector<int> &fault_main, std::vector<std::vector<int>> &fault_nodes_org,std::vector<int> &BIE_top, std::vector<int> &BIE_bot, std::vector<int> &AB_left, std::vector<int> &AB_right)
{
    std::string s;
    std::ifstream _in(fname);
    int mm = 0;
    while (true)
    {
        std::getline(_in,s);
        // Convert s to uppercase
        std::string upper(s);
        std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
        // 0.) Look for the "*Part" Section
//        if (upper.find("*PART")== static_cast<std::string::size_type>(0))
//        {
//            std::cout<<"Find *PART"<<std::endl;
//        }
        // 1.) Loo for the "*Nodes" section
        if (upper.find("*NODE")== static_cast<std::string::size_type>(0))
        {
            //std::string nset_name = s ;
            //std::cout<<s<<std::endl;
            // Temperatry variables for parsing lines of text
            char c ;
            std::string line;
            int i=0;
            while (_in.peek()!='*'&&_in.peek()!=EOF)
            {
                // Read an entire line which corresponds to a single points's id and (x,y) value
                std::getline(_in,line);
                // Revomie all whitesspaces characters from the line
                line.erase(std::remove_if(line.begin(),line.end(),isspace),line.end());
                // Make a stream out of the modified line so we can stream values from it in the usaly way
                std::stringstream ss(line);
                int abaqus_node_id = 0;
                double x = 0 , y = 0;
                ss >> abaqus_node_id >> c >> x >>c >> y ;
                Nodes.conservativeResize(i+1, 2);
                Nodes.row(i)<< x , y ;
                i++;
            }
            //  std::cout<< Nodes<< "\n"<< Nodes.rows()<<std::endl;
        }
        else if (upper.find("*ELEMENT,")==static_cast<std::string::size_type>(0))
        {
            //std::string elset_name = s;
            //std::cout<<s <<std::endl;
            char c ;
            std::string line;
            int i = 0;
            while (_in.peek()!='*'&&_in.peek()!=EOF)
            {
                // Read an entire line which corresponds to a single element's id and connectivity value (Q4_only)
                std::getline(_in,line);
                // Revomie all whitesspaces characters from the line
                line.erase(std::remove_if(line.begin(),line.end(),isspace),line.end());
                // Make a stream out of the modified line so we can stream values from it in the usaly way
                std::stringstream ss(line);
                int abaqus_el_id = 0;
                int node1 =0 , node2 =0 , node3 = 0 , node4 =0;
                ss >> abaqus_el_id >> c >> node1 >> c >>node2 >> c >> node3 >> c >> node4 ;
                // add -1 here becasue we are using a zero node numbering zero is the first node
                Element.conservativeResize(i+1, 4);
                Element.row(i)<<node1-1,node2-1,node3-1,node4-1;
                i++;
            }
            //   std::cout<<Element<<"\n" << Element.rows()<<std::endl;
        }
        else if (upper.find("*NSET, NSET=FAULT")==static_cast<std::string::size_type>(0))
        {
            std::string nset_name = s ;
            //std::cout<<s <<std::endl;
            //std::vector<int> fault_nodes;
            while (_in.peek()!='*'&&_in.peek()!=EOF)
            {
                // Read entire comma-seperated line into a string
                std::string csv_line;
                std::getline(_in, csv_line);
                // On that line use std::getline again to parse each comma-separated entry
                std::string cell;
                std::stringstream line_stream(csv_line);
                while (std::getline(line_stream,cell,','))
                {
                    char * endptr;
                    int id = static_cast<int>(std::strtol(cell.c_str(),&endptr,/*base=*/10));
                    // Note that lists of comma-spearated values in abaqus also
                    // "end" with a comma, so the last call to getline on a givenline will
                    // get an empty string, which we must delelet
                    if (id !=0 || cell.c_str() !=endptr)
                    {
                        // 'cell' is now a string with an integer id in it
                        // add -1 here becasue we are using a zero node numbering zero is the first node
                        fault_nodes_org[mm].push_back(id-1);
                    }
                }
            }
            mm = mm+1;
        }
        else if (upper.find("*NSET, NSET=MAIN_FAULT")==static_cast<std::string::size_type>(0))
        {
            std::string nset_name = s ;
            //std::cout<<s <<std::endl;
            //std::vector<int> fault_nodes;
            while (_in.peek()!='*'&&_in.peek()!=EOF)
            {
                // Read entire comma-seperated line into a string
                std::string csv_line;
                std::getline(_in, csv_line);
                // On that line use std::getline again to parse each comma-separated entry
                std::string cell;
                std::stringstream line_stream(csv_line);
                while (std::getline(line_stream,cell,','))
                {
                    char * endptr;
                    int id = static_cast<int>(std::strtol(cell.c_str(),&endptr,/*base=*/10));
                    // Note that lists of comma-spearated values in abaqus also
                    // "end" with a comma, so the last call to getline on a givenline will
                    // get an empty string, which we must delelet
                    if (id !=0 || cell.c_str() !=endptr)
                    {
                        // 'cell' is now a string with an integer id in it
                        // add -1 here becasue we are using a zero node numbering zero is the first node
                        fault_main.push_back(id-1);
                    }
                }
            }
        }
        else if (upper.find("*NSET, NSET=BIE_TOP")==static_cast<std::string::size_type>(0))
        {
            std::string nset_name = s ;
            //std::cout<<s <<std::endl;
            //std::vector<int> fault_nodes;
            while (_in.peek()!='*'&&_in.peek()!=EOF)
            {
                // Read entire comma-seperated line into a string
                std::string csv_line;
                std::getline(_in, csv_line);
                // On that line use std::getline again to parse each comma-separated entry
                std::string cell;
                std::stringstream line_stream(csv_line);
                while (std::getline(line_stream,cell,','))
                {
                    char * endptr;
                    int id = static_cast<int>(std::strtol(cell.c_str(),&endptr,/*base=*/10));
                    // Note that lists of comma-spearated values in abaqus also
                    // "end" with a comma, so the last call to getline on a givenline will
                    // get an empty string, which we must delelet
                    if (id !=0 || cell.c_str() !=endptr)
                    {
                        // 'cell' is now a string with an integer id in it
                        // add -1 here becasue we are using a zero node numbering zero is the first node
                        BIE_top.push_back(id-1);
                    }
                }
            }
        }
        else if (upper.find("*NSET, NSET=BIE_BOT")==static_cast<std::string::size_type>(0))
        {
            std::string nset_name = s ;
            //std::cout<<s <<std::endl;
            //std::vector<int> fault_nodes;
            while (_in.peek()!='*'&&_in.peek()!=EOF)
            {
                // Read entire comma-seperated line into a string
                std::string csv_line;
                std::getline(_in, csv_line);
                // On that line use std::getline again to parse each comma-separated entry
                std::string cell;
                std::stringstream line_stream(csv_line);
                while (std::getline(line_stream,cell,','))
                {
                    char * endptr;
                    int id = static_cast<int>(std::strtol(cell.c_str(),&endptr,/*base=*/10));
                    // Note that lists of comma-spearated values in abaqus also
                    // "end" with a comma, so the last call to getline on a givenline will
                    // get an empty string, which we must delelet
                    if (id !=0 || cell.c_str() !=endptr)
                    {
                        // 'cell' is now a string with an integer id in it
                        // add -1 here becasue we are using a zero node numbering zero is the first node
                        BIE_bot.push_back(id-1);
                    }
                }
            }
        }
        // Finding the Left and right Absorbing nodes
        else if (upper.find("*NSET, NSET=AB_LEFT")==static_cast<std::string::size_type>(0))
        {
            std::string nset_name = s ;
            while (_in.peek()!='*'&&_in.peek()!=EOF)
            {
                // Read entire comma-seperated line into a string
                std::string csv_line;
                std::getline(_in, csv_line);
                // On that line use std::getline again to parse each comma-separated entry
                std::string cell;
                std::stringstream line_stream(csv_line);
                while (std::getline(line_stream,cell,','))
                {
                    char * endptr;
                    int id = static_cast<int>(std::strtol(cell.c_str(),&endptr,/*base=*/10));
                    // Note that lists of comma-spearated values in abaqus also
                    // "end" with a comma, so the last call to getline on a givenline will
                    // get an empty string, which we must delelet
                    if (id !=0 || cell.c_str() !=endptr)
                    {
                        // 'cell' is now a string with an integer id in it
                        // add -1 here becasue we are using a zero node numbering zero is the first node
                        AB_left.push_back(id-1);
                    }
                }
            }
        }
        else if (upper.find("*NSET, NSET=AB_RIGHT")==static_cast<std::string::size_type>(0))
        {
            std::string nset_name = s ;
            while (_in.peek()!='*'&&_in.peek()!=EOF)
            {
                // Read entire comma-seperated line into a string
                std::string csv_line;
                std::getline(_in, csv_line);
                // On that line use std::getline again to parse each comma-separated entry
                std::string cell;
                std::stringstream line_stream(csv_line);
                while (std::getline(line_stream,cell,','))
                {
                    char * endptr;
                    int id = static_cast<int>(std::strtol(cell.c_str(),&endptr,/*base=*/10));
                    // Note that lists of comma-spearated values in abaqus also
                    // "end" with a comma, so the last call to getline on a givenline will
                    // get an empty string, which we must delelet
                    if (id !=0 || cell.c_str() !=endptr)
                    {
                        // 'cell' is now a string with an integer id in it
                        // add -1 here becasue we are using a zero node numbering zero is the first node
                        AB_right.push_back(id-1);
                    }
                }
            }
        }




        if (_in.eof())
            break;
    }
}
