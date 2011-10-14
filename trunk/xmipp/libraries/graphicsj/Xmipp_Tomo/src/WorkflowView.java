/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

/**
 * Custom display for Workflows (trees) of UserActions
 * 
 * Alternatives (for directed graphs): JPowerGraph, Java Universal Network/Graph Framework
 * 
 * Notes: it would be more appealing to display each UserAction in a single line
 * using hyperlinks, but the problem is that JTree by default only renders
 * the components of the tree nodes. For example, if you return a JButton in
 * the getTreeCellRendererComponent method of a custom TreeCellRenderer,
 * the tree displays the rendered image of the button. The button is not
 * really there, so it does not react to events.
 * 
 * Tables have the same problem. There are hacks for simulating buttons
 * (or any other Component like Combolists) with custom CellEditors, but we prefer to use a different approach:
 * display UserAction details as child nodes of the action
 * 
 * At the bottom of this file there is an example of customization in the
 * line of the appealing view
 * 
 */
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.TreeModelEvent;
import javax.swing.event.TreeModelListener;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;
import javax.swing.tree.*;

// TODO: change default icons - actions are neither folders nor files
// TODO: change details text field appeareance so it does not promote changing its contents,or
// display the parameters with a different look - maybe a table?
// TODO: simplify user's life: the workflow is autosaved every time a step is introduced or changed, and is autoloaded
// (or created if there is none yet) from the current dir. The user does not need to worry about loading/saving
// workflows, and hence does not need import/export neither
public class WorkflowView extends JPanel implements TreeSelectionListener, TreeModelListener{
	public static int PREFERED_WIDTH = 350;
	private static Dimension BUTTONS_PANEL_PREF_SIZE = new Dimension(PREFERED_WIDTH,150),
		TREE_PANEL_PREF_SIZE= new Dimension(PREFERED_WIDTH,350);
	private static String DETAILS_LABEL="Details", COMMENTS_LABEL="Comments";

	// @see valueChanged
	public static String ACTION_SELECTED="ActionSelected";

	// the controller must have access to the model. By default, this class acts as
	// a controller
	private Workflow model;
	
	private JTree tree;
	private Form formPanel;
	/**
	 * Create the panel with all its subcomponents
	 * @param root
	 * @param controller
	 */
	public WorkflowView(Workflow newModel,Object controller) {
		super();
		//super("WorkflowView");

		if(controller == null)
			controller = this;
		
        setBorder(BorderFactory.createCompoundBorder(
        		BorderFactory.createTitledBorder("Workflow"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
	
        model=newModel;
        model.addTreeModelListener(this);
		
		// use a layout that will stretch tree to panel size
		setLayout(new BorderLayout());

		// Create tree
		tree = new JTree(model);
		tree.addTreeSelectionListener(this);
		
		// Set line style
		tree.putClientProperty("JTree.lineStyle", "Angled");
		tree.setRowHeight(0);
		tree.setExpandsSelectedPaths(true);
		tree.setScrollsOnExpand(true);
		
		// tree.setSelectionPath(new TreePath(model.getRoot()));
		
		// Put tree in a scrollable pane
		JScrollPane sp = new JScrollPane(tree);
		sp.setPreferredSize(TREE_PANEL_PREF_SIZE);

		add(sp, BorderLayout.CENTER);
		formPanel=new Form("form",2);
		add(formPanel,BorderLayout.SOUTH);
		Action a = new Action(controller, new Command("workflow.print","Print Workflow","printWorkflow",true,null));
		// TODO: discard should delete the selected node and its children, with all the corresponding files
		Action b = new Action(controller, new Command("workflow.discardop","Discard Operation","discardOperation",true,null));
		// TODO: workflow autoload (from ./.project) then remove this button
		// TODO: workflow autosave, then remove the save button
		formPanel.addButton(a);
		formPanel.addButton(b);
		formPanel.addStringField(DETAILS_LABEL, null,null);
		formPanel.addStringField(COMMENTS_LABEL, null,getCommentsDocument());
		
		formPanel.setPreferredSize(BUTTONS_PANEL_PREF_SIZE);
	}

/**
 * @deprecated
 * @param action
 */
	public void newOperation(UserAction action){
		if(getSelectedNode() != null){
			DefaultMutableTreeNode n1=model.addUserAction(getSelectedNode(), action);
			tree.scrollPathToVisible(new TreePath(n1.getPath()));
		}
	}

	public void discardOperation(){
	    if (getSelectedNode() != null) {
            DefaultMutableTreeNode currentNode = getSelectedNode();
            MutableTreeNode parent = (MutableTreeNode)(currentNode.getParent());
            if (parent != null) {
            	model.removeNodeFromParent(currentNode);
            }
        } 
	}

	public void expandAll(){
		for (int i = 0; i < tree.getRowCount(); i++) 
	         tree.expandRow(i);
	}
	
	/**
	 * 
	 * @return the node selected in the tree (which may not be the same as in the model, until the model is updated to reflect it)
	 */
	public DefaultMutableTreeNode getSelectedNode(){
		return (DefaultMutableTreeNode) tree.getLastSelectedPathComponent();
	}
	
	/**
	 * @return may be null (watch out)
	 */
	public UserAction getSelectedUserAction(){
		if(getSelectedNode()==null)
			return null;
		return (UserAction) (getSelectedNode().getUserObject());
	}
	
	public void loadWorkflow(){
		Logger.debug("loadWorkflow");
		
	}
	
	public void printWorkflow(){
		Logger.debug(model.toString());
	}

	public void saveWorkflow(){
		Logger.debug("saveWorkflow");
		expandAll();
		TreeNode node= model.getSelectedNode();
		TreePath lastInsertedPath = new TreePath(model.getPathToRoot(node));
		Logger.debug(lastInsertedPath.toString());
		tree.setSelectionPath(lastInsertedPath);
		tree.scrollPathToVisible(lastInsertedPath);
		
	}
	
	private void setDetails(String details){
		if(details != null && formPanel != null)
			formPanel.setText(DETAILS_LABEL, details);
	}
	
	private void setCommentsDocument(Document doc){
		if(doc != null && formPanel != null)
			formPanel.setDocument(COMMENTS_LABEL, doc);
	}
	
	private Document getCommentsDocument(){
		UserAction ua=getSelectedUserAction();
		if(ua==null)
			return new PlainDocument();
		return ua.getCommentsDocument();
	}

	public void valueChanged(TreeSelectionEvent e){
		UserAction ua=getSelectedUserAction();
		if(ua != null){
			setDetails(ua.getCommandDetails());
			setCommentsDocument(ua.getCommentsDocument());
			model.setSelectedNode(getSelectedNode());
			// notify all listeners that the user selected an action in the workflow
			firePropertyChange(ACTION_SELECTED, false, true);
		}
	}
	
	public static void main(String args[]) {
		int windowWidth = 400, windowHeight = 300;
		JFrame window=new JFrame();
		Container content = window.getContentPane();
		Workflow workflow=Workflow.getTestWorkflow();

		WorkflowView wv = new WorkflowView(workflow,null);
		content.add(wv, BorderLayout.CENTER);
		window.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {System.exit(0);}
		});

		window.setSize(windowWidth, windowHeight);
		window.setVisible(true);
	}

	@Override
	public void treeNodesChanged(TreeModelEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void treeNodesInserted(TreeModelEvent e) {
		DefaultMutableTreeNode parent=(DefaultMutableTreeNode)(e.getTreePath().getLastPathComponent());
		TreeNode lastInserted = parent.getLastChild();
		TreePath lastInsertedPath = model.getNodePath(lastInserted);
		tree.setSelectionPath(lastInsertedPath);
		tree.scrollPathToVisible(lastInsertedPath);
	}

	@Override
	public void treeNodesRemoved(TreeModelEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void treeStructureChanged(TreeModelEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	
	
}

/*
 * Example of customization
 */
/* public class SwingLink extends JEditorPane implements HyperlinkListener{
private URI uri;

public SwingLink(){
	super("text/html","");
	setEditable(false);
    setOpaque (true);
    setPreferredSize(new Dimension(100,20));

	addHyperlinkListener(this);
}

public SwingLink(String text, URI uri){
    super();
    setup(text,uri);
}

public SwingLink(String text, String uri){
    super();
    URI oURI;
    try {
        oURI = new URI(uri);
    } catch (URISyntaxException e) {
        // converts to runtime exception for ease of use
        // if you cannot be sure at compile time that your
        // uri is valid, construct your uri manually and
        // use the other constructor.
        throw new RuntimeException(e);
    }
    setup(text,oURI);
}

public void hyperlinkUpdate(HyperlinkEvent event) {
    if (event.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
    	open(uri);
    }
}

public void setup(String t, URI u){
    uri = u;
    setText(t);
    setToolTipText(uri.toString());
}

@Override
public void setText(String text){
    setText(text,true);
}

public void setText(String text, boolean ul){
    String link = ul ? "<u>"+text+"</u>" : text;
    super.setText("<html>"+link+"</html>");
}

private void open(URI uri) {
	setText("Open" + getText());
}
}

private class CellRenderer implements TreeCellRenderer {

private SwingLink renderer;

CellRenderer () {
renderer = new SwingLink();
// setVerticalTextPosition(SwingConstants.TOP);
}

public Component  getTreeCellRendererComponent(JTree tree,
Object value, boolean selected, boolean expanded,
boolean leaf, int row, boolean hasFocus) {

// Change background color based on selected state
Color background = (selected ? Color.lightGray : Color.white);
renderer.setBackground(background);

UserAction action = (UserAction) ((DefaultMutableTreeNode) value).getUserObject();

String actionName = action.getName(), actionText="";
if(actionName != null)
  actionText = actionText + actionName + ".";

String actionCommand = action.getCommand();
if(actionCommand != null){
  actionText = actionText + actionCommand;
  String actionParameters = action.getParameters();
  if(actionParameters!= null)
	  actionText = actionText + " [Parameters]";
  if(action.getPlugin()!= null)
	  actionText = actionText + " [Options]";// getPlugin().getOptions();

  actionText = actionText + " :: " + action.getProgress();
}	

renderer.setText(actionText);

return renderer;
}
}*/