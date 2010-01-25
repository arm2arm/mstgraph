#include "program_settings.h"

void programm_settings::load(const std::string &filename) {
    using boost::property_tree::ptree;
    ptree pt;
    std::string ext(filename.end() - 3, filename.end());
    if (ext == "ini") {
        read_ini(filename, pt);
        m_test_type = pt.get<int>("OPTIONS.test_flag");
        if (m_test_type == 1)
            m_snap_file = pt.get<std::string > ("OPTIONS.test_file");
        else
            m_snap_file = pt.get<std::string > ("OPTIONS.file");

        m_path = pt.get<std::string > ("OPTIONS.path");
        m_HOP_file = pt.get<std::string > ("OPTIONS.HOP_file");
    } else {

        read_xml(filename, pt);
        int test_flag = pt.get("OPTIONS.test_flag", 0);

        BOOST_FOREACH(ptree::value_type &v,
                pt.get_child("OPTIONS.snaps"))
                m_snapvec.push_back(v.second.data());

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.models")) {
            m_modelvec.push_back(TPathMask(m.second.data(), "", ""));
        }
        if (m_modelvec.size() > 0) {
            int imodel = 0;

            BOOST_FOREACH(ptree::value_type &p,
                    pt.get_child("OPTIONS.paths")) {
                m_modelvec[imodel++].path = p.second.data();
            }
            imodel = 0;

            BOOST_FOREACH(ptree::value_type &p,
                    pt.get_child("OPTIONS.masks")) {
                m_modelvec[imodel++].mask = p.second.data();
            }

        }// end of modelvec.size()


        m_AnnTune = pt.get("OPTIONS.EST.ANNtuning", 12);
        m_HOPname.path = pt.get("OPTIONS.HOP.path", "");
        m_HOPname.mask = pt.get("OPTIONS.HOP.file", "");
        if (test_flag > 1) {
            m_TESTname.path = pt.get("OPTIONS.test.path", "");
            m_TESTname.mask = pt.get("OPTIONS.test.file", "");
        }

        int iest = 0;

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.EST.ndims")) {
            m_estvec.push_back(TEstData(m.second.get_value<int>()));
            std::cout << m_estvec[iest++].dim << std::endl;
        }

        iest = 0;

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.EST.ngbs")) {
            m_estvec[iest++].ngb = m.second.get_value<int>();
        }
        iest = 0;

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.EST.flags")) {
            m_estvec[iest++].flag = m.second.get_value<int>();
        }
        iest = 0;

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.EST.smooth_flags")) {
            m_estvec[iest++].smooth_flag = m.second.get_value<int>();
        }
        iest = 0;

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.EST.nsmooths")) {
            m_estvec[iest++].nsmooth = m.second.get_value<int>();
        }

        iest = 0;

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.EST.klinks")) {
            m_estvec[iest++].klink = m.second.get_value<int>();
        }
        iest = 0;

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.EST.files.rhofiles")) {
            m_estvec[iest++].rhofile = m_HOPname.path + m.second.data();
        }

        iest = 0;

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.EST.files.hsmlfiles")) {
            m_estvec[iest++].hsmlfile = m_HOPname.path + m.second.data();
        }
        iest = 0;

        BOOST_FOREACH(ptree::value_type &m,
                pt.get_child("OPTIONS.EST.files.ngbfiles")) {
            m_estvec[iest++].ngbfile = m_HOPname.path + m.second.data();
        }




    }
    print();
}

struct programm_settings pset;
