//===----- CFLSolver.cpp -- Context-free language reachability solver--------------//
//
//                     SVF: Static Value-Flow Analysis
//
// Copyright (C) <2013->  <Yulei Sui>
//

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//===----------------------------------------------------------------------===//

/*
 * CFLSolver.cpp
 *
 *  Created on: March 5, 2022
 *      Author: Yulei Sui
 */

#include "CFL/CFLSolver.h"
using namespace SVF;

double CFLSolver::numOfChecks = 0;

void CFLSolver::initialize()
{
    for(auto it = graph->begin(); it!= graph->end(); it++)
    {
        for(const CFLEdge* edge : (*it).second->getOutEdges())
        {
            pushIntoWorklist(edge);
        }
    }

    /// Foreach production X -> epsilon
    ///     add X(i,i) if not exist to E and to worklist
    for(const Production& prod : grammar->getEpsilonProds())
    {
        for(auto it = graph->begin(); it!= graph->end(); it++)
        {
            Symbol X = grammar->getLHSSymbol(prod);
            CFLNode* i = (*it).second;
            if(const CFLEdge* edge = graph->addCFLEdge(i, i, X))
            {
                pushIntoWorklist(edge);
            }
        }
    }
}

void CFLSolver::processCFLEdge(const CFLEdge* Y_edge)
{
    CFLNode* i = Y_edge->getSrcNode();
    CFLNode* j = Y_edge->getDstNode();

    /// For each production X -> Y
    ///     add X(i,j) if not exist to E and to worklist
    Symbol Y = Y_edge->getEdgeKind();
    if (grammar->hasProdsFromSingleRHS(Y))
        for(const Production& prod : grammar->getProdsFromSingleRHS(Y))
        {
            Symbol X = grammar->getLHSSymbol(prod);
            numOfChecks++;
            if(const CFLEdge* newEdge = graph->addCFLEdge(i, j, X))
            {
                pushIntoWorklist(newEdge);
            }
        }

    /// For each production X -> Y Z
    /// Foreach outgoing edge Z(j,k) from node j do
    ///     add X(i,k) if not exist to E and to worklist
    if (grammar->hasProdsFromFirstRHS(Y))
        for(const Production& prod : grammar->getProdsFromFirstRHS(Y))
        {
            Symbol X = grammar->getLHSSymbol(prod);
            for(const CFLEdge* Z_edge : j->getOutEdgeWithTy(grammar->getSecondRHSSymbol(prod)))
            {
                CFLNode* k = Z_edge->getDstNode();
                numOfChecks++;
                if(const CFLEdge* newEdge = graph->addCFLEdge(i, k, X))
                {
                    pushIntoWorklist(newEdge);
                }
            }
        }

    /// For each production X -> Z Y
    /// Foreach incoming edge Z(k,i) to node i do
    ///     add X(k,j) if not exist to E and to worklist
    if(grammar->hasProdsFromSecondRHS(Y))
        for(const Production& prod : grammar->getProdsFromSecondRHS(Y))
        {
            Symbol X = grammar->getLHSSymbol(prod);
            for(const CFLEdge* Z_edge : i->getInEdgeWithTy(grammar->getFirstRHSSymbol(prod)))
            {
                CFLNode* k = Z_edge->getSrcNode();
                numOfChecks++;
                if(const CFLEdge* newEdge = graph->addCFLEdge(k, j, X))
                {
                    pushIntoWorklist(newEdge);
                }
            }
        }
}


void CFLSolver::solve()
{
    /// initial worklist
    initialize();

    while(!isWorklistEmpty())
    {
        /// Select and remove an edge Y(i,j) from worklist
        const CFLEdge* Y_edge = popFromWorklist();
        processCFLEdge(Y_edge);
    }
}

void POCRSolver::buildCFLData()
{
    for (CFLEdge* edge: graph->getCFLEdges())
        addEdge(edge->getSrcID(), edge->getDstID(), edge->getEdgeKind());
}

void POCRSolver::processCFLEdge(const CFLEdge* Y_edge)
{
    CFLNode* i = Y_edge->getSrcNode();
    CFLNode* j = Y_edge->getDstNode();

    /// For each production X -> Y
    ///     add X(i,j) if not exist to E and to worklist
    Symbol Y = Y_edge->getEdgeKind();
    if (grammar->hasProdsFromSingleRHS(Y))
        for(const Production& prod : grammar->getProdsFromSingleRHS(Y))
        {
            Symbol X = grammar->getLHSSymbol(prod);
            numOfChecks++;
            if (addEdge(i->getId(), j->getId(), X))
            {
                const CFLEdge* newEdge = graph->addCFLEdge(Y_edge->getSrcNode(), Y_edge->getDstNode(), X);
                pushIntoWorklist(newEdge);
            }

        }

    /// For each production X -> Y Z
    /// Foreach outgoing edge Z(j,k) from node j do
    ///     add X(i,k) if not exist to E and to worklist
    if (grammar->hasProdsFromFirstRHS(Y))
        for(const Production& prod : grammar->getProdsFromFirstRHS(Y))
        {
            Symbol X = grammar->getLHSSymbol(prod);
            NodeBS diffDsts = addEdges(i->getId(), getSuccMap(j->getId())[grammar->getSecondRHSSymbol(prod)], X);
            numOfChecks += getSuccMap(j->getId())[grammar->getSecondRHSSymbol(prod)].count();
            for (NodeID diffDst: diffDsts)
            {
                const CFLEdge* newEdge = graph->addCFLEdge(i, graph->getGNode(diffDst), X);
                pushIntoWorklist(newEdge);
            }
        }

    /// For each production X -> Z Y
    /// Foreach incoming edge Z(k,i) to node i do
    ///     add X(k,j) if not exist to E and to worklist
    if(grammar->hasProdsFromSecondRHS(Y))
        for(const Production& prod : grammar->getProdsFromSecondRHS(Y))
        {
            Symbol X = grammar->getLHSSymbol(prod);
            NodeBS diffSrcs = addEdges(getPredMap(i->getId())[grammar->getFirstRHSSymbol(prod)], j->getId(), X);
            numOfChecks += getPredMap(i->getId())[grammar->getFirstRHSSymbol(prod)].count();
            for (NodeID diffSrc: diffSrcs)
            {
                const CFLEdge* newEdge = graph->addCFLEdge(graph->getGNode(diffSrc), j, X);
                pushIntoWorklist(newEdge);
            }
        }
}

void POCRSolver::initialize()
{
    for(auto edge : graph->getCFLEdges())
    {
        pushIntoWorklist(edge);
    }

    /// Foreach production X -> epsilon
    ///     add X(i,i) if not exist to E and to worklist
    for(const Production& prod : grammar->getEpsilonProds())
    {
        for(auto IDMap : getSuccMap())
        {
            Symbol X = grammar->getLHSSymbol(prod);
            if (addEdge(IDMap.first, IDMap.first, X))
            {
                CFLNode* i = graph->getGNode(IDMap.first);
                const CFLEdge* newEdge = graph->addCFLEdge(i, i, X);
                pushIntoWorklist(newEdge);
            }
        }
    }
}

void POCRHybridSolver::processCFLEdge(const CFLEdge* Y_edge)
{
    CFLNode* i = Y_edge->getSrcNode();
    CFLNode* j = Y_edge->getDstNode();

    /// For each production X -> Y
    ///     add X(i,j) if not exist to E and to worklist
    Symbol Y = Y_edge->getEdgeKind();
    if (grammar->hasProdsFromSingleRHS(Y))
        for(const Production& prod : grammar->getProdsFromSingleRHS(Y))
        {
            Symbol X = grammar->getLHSSymbol(prod);
            numOfChecks++;
            if (addEdge(i->getId(), j->getId(), X))
            {
                const CFLEdge* newEdge = graph->addCFLEdge(Y_edge->getSrcNode(), Y_edge->getDstNode(), X);
                pushIntoWorklist(newEdge);
            }

        }

    /// For each production X -> Y Z
    /// Foreach outgoing edge Z(j,k) from node j do
    ///     add X(i,k) if not exist to E and to worklist
    if (grammar->hasProdsFromFirstRHS(Y))
        for(const Production& prod : grammar->getProdsFromFirstRHS(Y))
        {
            if ((grammar->getLHSSymbol(prod) == grammar->strToSymbol("F")) && (Y == grammar->strToSymbol("F")) && (grammar->getSecondRHSSymbol(prod) == grammar->strToSymbol("F")))
            {
                addArc(i->getId(), j->getId());
            }
            else
            {
                Symbol X = grammar->getLHSSymbol(prod);
                NodeBS diffDsts = addEdges(i->getId(), getSuccMap(j->getId())[grammar->getSecondRHSSymbol(prod)], X);
                numOfChecks += getSuccMap(j->getId())[grammar->getSecondRHSSymbol(prod)].count();
                for (NodeID diffDst: diffDsts)
                {
                    const CFLEdge* newEdge = graph->addCFLEdge(i, graph->getGNode(diffDst), X);
                    pushIntoWorklist(newEdge);
                }
            }
        }

    /// For each production X -> Z Y
    /// Foreach incoming edge Z(k,i) to node i do
    ///     add X(k,j) if not exist to E and to worklist
    if(grammar->hasProdsFromSecondRHS(Y))
        for(const Production& prod : grammar->getProdsFromSecondRHS(Y))
        {
            if ((grammar->getLHSSymbol(prod) == grammar->strToSymbol("F")) && (Y == grammar->strToSymbol("F")) && (grammar->getFirstRHSSymbol(prod) == grammar->strToSymbol("F")))
            {
                addArc(i->getId(), j->getId());
            }
            else
            {
                Symbol X = grammar->getLHSSymbol(prod);
                NodeBS diffSrcs = addEdges(getPredMap(i->getId())[grammar->getFirstRHSSymbol(prod)], j->getId(), X);
                numOfChecks += getPredMap(i->getId())[grammar->getFirstRHSSymbol(prod)].count();
                for (NodeID diffSrc: diffSrcs)
                {
                    const CFLEdge* newEdge = graph->addCFLEdge(graph->getGNode(diffSrc), j, X);
                    pushIntoWorklist(newEdge);
                }
            }
        }
}

void POCRHybridSolver::initialize()
{
    for(auto edge : graph->getCFLEdges())
    {
        pushIntoWorklist(edge);
    }

    // init hybrid dataset
    for (auto it = graph->begin(); it != graph->end(); ++it)
    {
        NodeID nId = it->first;
        addInd_h(nId, nId);
    }

    ///     add X(i,i) if not exist to E and to worklist
    for(const Production& prod : grammar->getEpsilonProds())
    {
        for(auto IDMap : getSuccMap())
        {
            Symbol X = grammar->getLHSSymbol(prod);
            if (addEdge(IDMap.first, IDMap.first, X))
            {
                CFLNode* i = graph->getGNode(IDMap.first);
                const CFLEdge* newEdge = graph->addCFLEdge(i, i, X);
                pushIntoWorklist(newEdge);
            }
        }
    }
}

void POCRHybridSolver::addArc(NodeID src, NodeID dst)
{
    if(hasEdge(src, dst, grammar->strToSymbol("F")))
        return;

    for (auto& iter: indMap[src])
    {
        meld(iter.first, getNode_h(iter.first, src), getNode_h(dst, dst));
    }
}


void POCRHybridSolver::meld(NodeID x, TreeNode* uNode, TreeNode* vNode)
{
    numOfChecks++;

    TreeNode* newVNode = addInd_h(x, vNode->id);
    if (!newVNode)
        return;

    insertEdge_h(uNode, newVNode);
    for (TreeNode* vChild: vNode->children)
    {
        meld_h(x, newVNode, vChild);
    }
}

bool POCRHybridSolver::hasInd_h(NodeID src, NodeID dst)
{
    auto it = indMap.find(dst);
    if (it == indMap.end())
        return false;
    return (it->second.find(src) != it->second.end());
}

POCRHybridSolver::TreeNode* POCRHybridSolver::addInd_h(NodeID src, NodeID dst)
{
    TreeNode* newNode = new TreeNode(dst);
    auto resIns = indMap[dst].insert(std::make_pair(src, newNode));
    if (resIns.second)
        return resIns.first->second;
    delete newNode;
    return nullptr;
}

void POCRHybridSolver::addArc_h(NodeID src, NodeID dst)
{
    if (!hasInd_h(src, dst))
    {
        for (auto iter: indMap[src])
        {
            meld_h(iter.first, getNode_h(iter.first, src), getNode_h(dst, dst));
        }
    }
}

void POCRHybridSolver::meld_h(NodeID x, TreeNode* uNode, TreeNode* vNode)
{
    TreeNode* newVNode = addInd_h(x, vNode->id);
    if (!newVNode)
        return;

    insertEdge_h(uNode, newVNode);
    for (TreeNode* vChild: vNode->children)
    {
        meld_h(x, newVNode, vChild);
    }
}


// На вход подаётся набор kind-ов: terminals/nonterminals.
// По ним строится упорядоченное множество самих терминалов и нетерминалов уже с учётом аттрибутов.
// Returns size of filled enumerated_symbols. 
/*uint64_t MatrixSolver::enumerate(Map<std::string, Kind> kinds, CFGrammar::SymbolMap<Label, uint32_t>& enumerated_symbols){
    uint32_t index = 0;
    Set<Kind> attributedKinds = grammar->getAttrSyms();
    for (auto str2kind: kinds){
        Label sym;
        Kind kind = str2kind.second;
        sym.kind = kind;
        // printf("Kind: %s\n",(grammar->kindToStr(kind)).c_str());
        for(auto [key, value]:grammar->getKindToAttrsMap()){
            printf("Kind: %s\n",(grammar->kindToStr(key)).c_str());
            for(auto i:value){
                printf(" attr: %d ", i);
            }
        }
        printf("grammar->getAttrSyms():\n");
        for(auto vkind: grammar->getAttrSyms()){
            printf("Kind: %s\n",(grammar->kindToStr(vkind)).c_str());
        }
        if (attributedKinds.find(kind) != attributedKinds.end()){
            for(auto attri: grammar->getKindToAttrsMap().at(kind)){
                sym.attribute = attri;
                enumerated_symbols[sym] = index++;
            }
        } else{
            enumerated_symbols[sym] = index++;
        }
    }
    return index;
}
*/

void MatrixSolver::enumerate(){

    //grammar->dump();
    uint32_t nonterm_total = 0;
    uint32_t term_total = 0;
    for(auto pro: grammar->getEpsilonProds()){
        if(enumerated_nonterminals.find(pro[0]) == enumerated_nonterminals.end()){
            enumerated_nonterminals[pro[0]] = nonterm_total++;
        }
    }

    for (auto sym2prods: grammar->getSingleRHSToProds())
    {
       for (auto pro: sym2prods.second) {
        if(enumerated_nonterminals.find(pro[0]) == enumerated_nonterminals.end()){
            enumerated_nonterminals[pro[0]] = nonterm_total++;
        }
        if(enumerated_terminals.find(pro[1]) == enumerated_terminals.end()){
            enumerated_terminals[pro[1]] = term_total++;
        }
       }
    }

    for (auto sym2prods: grammar->getFirstRHSToProds())
    {
        for (auto pro: sym2prods.second) {

            for (int i = 0; i<3; i++){
                if ( (grammar->getTerminals()).find(grammar->kindToStr(pro[i].kind)) != (grammar->getTerminals()).end() ){
                    if(enumerated_terminals.find(pro[i]) == enumerated_terminals.end()){
                        enumerated_terminals[pro[i]] = term_total++;
                    }
                } else {
                    if(enumerated_nonterminals.find(pro[i]) == enumerated_nonterminals.end()){
                        enumerated_nonterminals[pro[i]] = nonterm_total++;
                    }
                }
            }
        }
        
    }
    printf("Terminals:\n");
    for (auto [term, index]: enumerated_terminals){
        printf("Term %d is %s\n",index, (grammar->symToStrDump(term)).c_str());
    }
    printf("Nonterminals:\n");
    for (auto [nonterm, index]: enumerated_nonterminals){
        printf("Nonterm %d is %s\n",index, (grammar->symToStrDump(nonterm)).c_str());
    }
}

void MatrixSolver::graphSVF2LAGraph(GrB_Matrix *adj_matrices, int64_t totalTerm)
{
    typedef std::pair<std::vector<uint64_t>, std::vector<uint64_t>> VectorPair;
    std::vector<VectorPair> nodeID_pairs(totalTerm);
    uint64_t totalNode = graph->getTotalNodeNum(); // uint32_t -> uint64_t

    for (auto edge: graph->getCFLEdges()){
        /*Label sym;
        sym.kind = 0;
        sym.attribute = edge->getEdgeAttri();
        printf("EdgeFlag is %lld\n", edge->getEdgeKind());
        sym.variableAttribute = 0;
        printf("Symbol is %s\n", (grammar->symToStrDump(sym)).c_str());
        uint64_t term_index = enumerated_terminals.at(edge->getEdgeKind()); //FIX (abort)
        printf("Тут\n");
        nodeID_pairs[term_index].first.push_back(edge->getSrcNode()->getId());
        nodeID_pairs[term_index].second.push_back(edge->getDstNode()->getId());*/
        printf(" %s ", (grammar->symToStrDump(edge->getEdgeKind())).c_str());
    }
    exit(0);

    for (int i = 0; i < totalTerm; ++i){
        uint64_t* row_indexes = (nodeID_pairs[i]).first.data();
        uint64_t* column_indexes = (nodeID_pairs[i]).second.data();
        uint64_t nvals = (nodeID_pairs[i]).first.size();
        bool values[nvals];
        std::fill_n(values, nvals, true);
        GRB_TRY(GrB_Matrix_new(&adj_matrices[i], GrB_BOOL, totalNode, totalNode));
        GRB_TRY(GrB_Matrix_build_BOOL(adj_matrices[i], row_indexes, column_indexes, values, nvals, nullptr)); // uint32_t -> uint64_t в индексах
    }

}

void MatrixSolver::grammarSVF2LAGraph(LAGraph_rule_WCNF *rules)
{
    int index = 0;
    for(auto pro: grammar->getEpsilonProds()){
        rules[index].nonterm = enumerated_nonterminals.at(pro[0]); // uint32 -> int32
        rules[index].prod_A = -1;
        rules[index].prod_B = -1;
        rules[index].index = 0;
        index++;
    }

    for (auto sym2prods: grammar->getSingleRHSToProds())
    {
        for (auto pro: sym2prods.second)
        {
            rules[index].nonterm = enumerated_nonterminals.at(pro[0]);
            rules[index].prod_A = enumerated_terminals.at(pro[1]);
            rules[index].prod_B = -1;
            rules[index].index = 0;
            index++;
        }
    }

    for (auto sym2prods: grammar->getFirstRHSToProds())
    {
        for (auto pro: sym2prods.second)
        {
            rules[index].nonterm = enumerated_nonterminals.at(pro[0]);
            rules[index].prod_A = enumerated_nonterminals.at(pro[1]);
            rules[index].prod_B = enumerated_nonterminals.at(pro[2]);
            rules[index].index = 0;
            index++;
        }
    }
}

void MatrixSolver::graphLAGraph2SVF(GrB_Matrix *nonterm_matrices, int64_t totalNonterm)
{
    for(int i = 0; i < totalNonterm; i++){

        uint64_t nvals;
        GRB_TRY(GrB_Matrix_nvals(&nvals, nonterm_matrices[i]));
        uint64_t row_indexes[nvals];
        uint64_t column_indexes[nvals];
        bool values[nvals];
        GRB_TRY(GrB_Matrix_extractTuples_BOOL(row_indexes, column_indexes, values, &nvals, nonterm_matrices[i]));
        for(uint64_t j = 0; j < nvals; j++){
            printf("Здесь\n");
            graph->addCFLEdge(graph->getGNode(row_indexes[j]), graph->getGNode(column_indexes[j]), enumerated_nonterminals.at(i));
            printf("Тут\n");
        }
    }

}

void MatrixSolver::solve()
{
    enumerate();
    int64_t totalTerm = enumerated_terminals.size();
    int64_t totalNonterm = enumerated_nonterminals.size();
    int64_t totalRules = 0;
    for (auto sym2prods: grammar->getRawProductions()){
        totalRules += sym2prods.second.size();
    }

    //init LAGraph objects
    GrB_Matrix adj_matrices[totalTerm];
    GrB_Matrix new_nonterm_edges[totalNonterm];
    LAGraph_rule_WCNF rules[totalRules];
    char* msg = NULL;
    graphSVF2LAGraph(adj_matrices, totalTerm);
    grammarSVF2LAGraph(rules);
    LAGRAPH_TRY(LAGraph_CFL_reachability(new_nonterm_edges, adj_matrices, totalTerm, totalNonterm, rules, totalRules, msg));
    if (msg != nullptr){
        assert(false && msg);
        abort();
    }
    graphLAGraph2SVF(new_nonterm_edges, totalNonterm);

}
