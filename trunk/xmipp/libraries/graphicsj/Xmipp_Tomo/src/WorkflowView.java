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
import javax.swing.tree.*;

// TODO: test insert/delete nodes in runtime
// TODO: should show/expand all nodes by default (except the details nodes)
// TODO: new design: add a form in the bottom for details & comments
public class WorkflowView extends JPanel {
	// The initial width and height of the frame
	private static int WIDTH = 300, HEIGHT = 200;
	private static Dimension BUTTONS_PANEL_PREF_SIZE = new Dimension(250,100);

	// the controller must have access to the model
	private Workflow model;
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
		
		model=newModel;
		TreeNode root = model.getRoot();
		
		// use a layout that will stretch tree to panel size
		setLayout(new BorderLayout());

		// Create tree
		JTree tree = new JTree(root);

		// Set line style
		tree.putClientProperty("JTree.lineStyle", "Angled");
		tree.setRowHeight(0);
		tree.setScrollsOnExpand(true);
		/*tree.setAlignmentX(Component.TOP_ALIGNMENT);
      tree.setAlignmentY(Component.TOP_ALIGNMENT);*/


		// Put tree in a scrollable pane
		JScrollPane sp = new JScrollPane(tree);

		add(sp, BorderLayout.CENTER);
		JPanel formPanel=new JPanel();
		formPanel.setLayout(new GridBagLayout());
		add(formPanel,BorderLayout.SOUTH);
		Action a = new Action(controller, new Command("workflow.newop","New Operation","newOperation",true,null));
		Action b = new Action(controller, new Command("workflow.discardop","Discard Operation","discardOperation",true,null));      
		Action c = new Action(controller, new Command("workflow.load","Load Workflow","loadWorkflow",true,null));
		Action d = new Action(controller, new Command("workflow.save","Save Workflow","saveWorkflow",true,null));
		JButton newOp=new JButton(a), delOp=new JButton(b), load=new JButton(c), save=new JButton(d);
		formPanel.add(newOp);
		formPanel.add(delOp);
		formPanel.add(load);
		formPanel.add(save);
		formPanel.setPreferredSize(BUTTONS_PANEL_PREF_SIZE);
	}

	// TODO: should add the operation to the current node
	// TODO: the treepanel does not  show new node...
	public void newOperation(){
		Logger.debug("newOperation");
		UserAction a1=new UserAction(0,"Load2","Load2" , "g1ta2.spi");
		DefaultMutableTreeNode n1=model.addUserAction(model.getRoot(), a1);
	}

	public void discardOperation(){
		Logger.debug("discardOperation");
	}

	public void loadWorkflow(){
		Logger.debug("loadWorkflow");
		Logger.debug(model.toString());
	}

	public void saveWorkflow(){
		Logger.debug("saveWorkflow");
	}

	public static void main(String args[]) {
		JFrame window=new JFrame();
		Container content = window.getContentPane();
		Workflow workflow=Workflow.getTestWorkflow();

		WorkflowView wv = new WorkflowView(workflow,null);
		content.add(wv, BorderLayout.CENTER);
		window.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {System.exit(0);}
		});

		window.setSize(WIDTH, HEIGHT);
		window.setVisible(true);
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